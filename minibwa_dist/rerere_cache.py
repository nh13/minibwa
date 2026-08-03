"""Persist `.git/rr-cache` on an orphan branch so CI runners start warm.

rerere resolutions are the difference between a sync that is free and a sync
that needs a human. The cache is local to a clone and is never pushed, so a
fresh Actions runner would re-hit every conflict from scratch. Storing it as a
tree on an orphan branch makes it durable and reviewable.

Runbook -- removing a bad entry: `save()` only ever adds entries; nothing here
ever deletes one. If a rerere resolution is semantically wrong -- the merge
looks clean, rerere replays it with no conflict, but the byte-identity gate
(`minibwa_dist/gates.py`) then fails -- that bad entry is now permanent: it will be
staged and replayed on every subsequent sync, forever, and the gate will fail
every single time until a human intervenes. There is no code path that detects
or expires a bad entry on its own. The remedy is manual: check out the
`rerere-cache` branch, delete the offending `<hash>/` directory by hand (the
gate's failure message names the file; cross-reference against the sync run's
merge commit to find which feature's resolution is at fault), commit, and push
`rerere-cache`. The next sync then hits the conflict fresh instead of
replaying the bad resolution.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from minibwa_dist.gitutil import git, git_ok, stdout

CACHE_BRANCH = "rerere-cache"


class OrphanOverPublishedCache(RuntimeError):
    """`save()` was asked to build a first-run orphan on top of a published cache."""


def _refuse_orphan_over_published_cache(repo: Path, branch: str) -> None:
    """Refuse the one shape of `save()` that can destroy the cache.

    With no local `refs/heads/<branch>`, `save()` builds a fresh orphan holding
    only the entries this run happened to record. That is right on a genuine
    first run and catastrophic after a *failed* fetch: the caller then publishes
    that orphan over a branch carrying every resolution ever recorded, wiping
    them from the tip and from the history alike, with no local signal that
    anything was lost.

    Locally the two cases differ in exactly one way -- a remote-tracking ref. CI
    mirrors every one of the fork's heads into `refs/remotes/origin/*` before it
    fetches the cache branch, so `refs/remotes/*/<branch>` exists whenever the
    published branch does. Local branch missing while the tracking ref is
    present means the restore did not happen, not that there is nothing to
    restore. Refuse loudly rather than manufacture a replacement: an
    unrestored cache means this run also re-hit every conflict from scratch,
    which is a build worth failing.
    """
    tracked = stdout(repo, "for-each-ref", "--format=%(refname)", f"refs/remotes/*/{branch}")
    if tracked:
        raise OrphanOverPublishedCache(
            f"refusing to build a fresh orphan '{branch}': {tracked.splitlines()[0]} exists, "
            f"so the cache is published but was never fetched into refs/heads/{branch}. "
            "Fetch it and re-run; do not publish over it."
        )


def _rr_cache_dir(repo: Path) -> Path:
    """Absolute path to this repo's rerere cache directory.

    Must resolve to the *common* gitdir, not the per-worktree one: rerere
    always writes `rr-cache` under the common gitdir, even when invoked from
    a linked worktree (`--absolute-git-dir` would instead return the private
    `.git/worktrees/<name>` directory there, which never holds the cache).
    """
    return (
        Path(stdout(repo, "rev-parse", "--path-format=absolute", "--git-common-dir")) / "rr-cache"
    )


def restore(repo: Path, branch: str = CACHE_BRANCH) -> int:
    """Populate `.git/rr-cache` from `branch`. Returns the number of entries.

    Missing branch is not an error — a first run legitimately has no cache.
    """
    if not git_ok(repo, "rev-parse", "--verify", f"refs/heads/{branch}"):
        return 0

    cache = _rr_cache_dir(repo)
    cache.mkdir(parents=True, exist_ok=True)
    entries: list[Path] = []  # bound before the `with`, so a worktree failure
    with tempfile.TemporaryDirectory() as tmp:  # surfaces as itself, not NameError
        git(repo, "worktree", "add", "--detach", tmp, branch)
        try:
            source = Path(tmp)
            entries = [d for d in source.iterdir() if d.is_dir() and d.name != ".git"]
            for entry in entries:
                shutil.copytree(entry, cache / entry.name, dirs_exist_ok=True)
        finally:
            git(repo, "worktree", "remove", "--force", tmp, check=False)
    return len(entries)


def save(repo: Path, branch: str = CACHE_BRANCH) -> bool:
    """Commit `.git/rr-cache` to `branch`. Returns True if anything changed.

    Raises `OrphanOverPublishedCache` when the branch is missing locally but a
    remote-tracking ref says it is published -- see that guard for why.
    """
    cache = _rr_cache_dir(repo)
    entries = sorted(d for d in cache.iterdir() if d.is_dir()) if cache.is_dir() else []
    if not entries:
        return False

    has_branch = git_ok(repo, "rev-parse", "--verify", f"refs/heads/{branch}")
    if not has_branch:
        _refuse_orphan_over_published_cache(repo, branch)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        if has_branch:
            git(repo, "worktree", "add", str(work), branch)
        else:
            git(repo, "worktree", "add", "--detach", str(work))
        # Everything from here on runs against a registered worktree, so the
        # `finally` below must cover it too -- otherwise a failure partway
        # through orphan setup (or anything after) leaves a stale, prunable
        # worktree registration behind.
        try:
            if not has_branch:
                git(work, "checkout", "--orphan", branch)
                git(work, "rm", "-rf", "--ignore-unmatch", ".")
            for entry in entries:
                shutil.copytree(entry, work / entry.name, dirs_exist_ok=True)
            # Stage by explicit path: `entries` already names every cache
            # directory, so there is no need to reach for `git add --all`.
            git(work, "add", "--", *(entry.name for entry in entries))
            if not git(work, "diff", "--cached", "--quiet", check=False).returncode:
                return False
            git(
                work,
                "-c",
                "commit.gpgsign=false",
                "commit",
                "-m",
                f"chore(rerere): {len(entries)} resolution(s)",
            )
        finally:
            git(repo, "worktree", "remove", "--force", str(work), check=False)
    return True
