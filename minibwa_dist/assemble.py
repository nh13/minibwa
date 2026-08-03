"""Build `dist-next` by merging frozen feature branches onto a fresh base.

Merge, never rebase. Rebasing replays a feature's commits onto the new
upstream, which shifts the content that rerere keys on and defeats the cache.
Merging leaves each feature byte-stable, so a feature-vs-feature conflict
presents the same pre-image every sync and replays automatically. Only
feature-vs-upstream conflicts can be novel.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from minibwa_dist.gitutil import git, git_ok, stdout
from minibwa_dist.manifest import Feature, Manifest


class AssemblyError(RuntimeError):
    """Anything that stops an assembly. Callers report these; they never traceback."""


class RequiredFeatureConflict(AssemblyError):
    """A `required = true` feature would not merge. The build must not proceed."""


class TrainingStop(RequiredFeatureConflict):
    """`--stop-at` reached its feature and left the conflict for a human.

    A subclass, not a sibling: callers that only know `RequiredFeatureConflict`
    still catch it, and the two outcomes -- a successful training stop and a
    failed build -- are now told apart by type instead of by matching the
    exception's message text, which any edit to either message would break.
    """


class MissingBranch(AssemblyError):
    """A manifest branch does not exist. Branch custody violation — never a silent drop."""


class MergeFailed(AssemblyError):
    """`git merge` refused to start. Not a conflict, and not something to commit."""


class UnknownStopTarget(AssemblyError):
    """`--stop-at` named something the manifest does not contain."""


@dataclass(frozen=True)
class DroppedFeature:
    """An optional feature excluded from this build because it would not merge."""

    name: str
    conflicting_files: tuple[str, ...]


@dataclass(frozen=True)
class AssemblyResult:
    """Outcome of one assembly run."""

    base: str
    head: str
    merged: tuple[str, ...]
    dropped: tuple[DroppedFeature, ...]


def _conflicting_files(repo: Path) -> tuple[str, ...]:
    out = stdout(repo, "diff", "--name-only", "--diff-filter=U")
    return tuple(line for line in out.splitlines() if line)


def _current_ref(repo: Path) -> str:
    """The ref to return the worktree to when the assembly is over.

    `rev-parse --abbrev-ref HEAD` answers the literal string `HEAD` on a
    detached checkout, so restoring it is `git checkout HEAD` -- a no-op that
    silently leaves the caller parked on the assembled branch. Ask for the
    symbolic ref, which fails cleanly when detached, and fall back to the commit.
    """
    branch = git(repo, "symbolic-ref", "-q", "--short", "HEAD", check=False).stdout.strip()
    return branch or stdout(repo, "rev-parse", "HEAD")


def _try_merge(repo: Path, feature: Feature, ref: str) -> tuple[bool, tuple[str, ...]]:
    """Merge `ref` for `feature`. Returns (merged, conflicting_files).

    `commit.gpgsign=false` is forced: `git merge` honours the setting, and with
    1Password signing that means 14 biometric prompts per run, or a hard failure
    when the key is unavailable. Bot-built integration merges on `dist` are
    unsigned by design; the signing rule governs human commits.

    rerere.autoUpdate stages any resolution it recognises, so a merge that exits
    nonzero but leaves no unmerged paths *may* be one rerere already handled --
    but the same signature also covers a merge git refused to begin (a dirty
    worktree, or an untracked file the merge would overwrite). `MERGE_HEAD`
    tells the two apart positively: it exists only once the merge actually
    started. Inferring the rerere hit from the absence of conflicts instead used
    to run `git commit --no-edit` on a repository that was never merging, which
    failed with git's own message swallowed.

    All three rerere/conflict-style knobs are pinned here, on the merge
    invocation itself, rather than left to ambient config, so this function's
    behaviour does not depend on whoever's environment happens to run it.
    `rerere.enabled` defaults to false until `.git/rr-cache` exists, which is
    exactly the state a fresh clone is in, so without the pin the bootstrap
    training records nothing. Without `rerere.autoUpdate`, rerere resolves the
    content but leaves the path unmerged, so `_conflicting_files` reports a
    conflict and every rerere *hit* reads as a drop instead of a merge.

    `merge.conflictStyle` is a weaker case, and the comment here used to
    overstate it: measured on git 2.50.1, rerere normalises its pre-image, so a
    resolution trained under an ambient `zdiff3` replays fine against a plain
    `merge`-style merge and vice versa. What the pin does buy is what the human
    sees during `--stop-at` training: the conflict left in the worktree carries
    the same markers CI would have produced, rather than whichever style that
    maintainer happens to prefer.
    """
    merge = git(
        repo,
        "-c",
        "rerere.enabled=true",
        "-c",
        "rerere.autoUpdate=true",
        "-c",
        "merge.conflictStyle=merge",
        "-c",
        "commit.gpgsign=false",
        "merge",
        "--no-ff",
        "--no-edit",
        "-m",
        f"merge {feature.name} ({feature.branch})",
        ref,
        check=False,
    )
    if merge.returncode == 0:
        return True, ()

    conflicts = _conflicting_files(repo)
    if not conflicts:
        if not git_ok(repo, "rev-parse", "-q", "--verify", "MERGE_HEAD"):
            raise MergeFailed(
                f"merging '{ref}' for feature '{feature.name}' never started: "
                f"{merge.stderr.strip() or merge.stdout.strip()}"
            )
        git(repo, "-c", "commit.gpgsign=false", "commit", "--no-edit")
        return True, ()

    return False, conflicts


def assemble(
    repo: Path,
    manifest: Manifest,
    base_ref: str,
    out_branch: str,
    remote: str,
    stop_at: str | None = None,
) -> AssemblyResult:
    """Assemble `out_branch` from `base_ref` plus every manifest feature.

    Features resolve as `{remote}/{branch}` -- `nh13/...` locally, `origin/...`
    in CI. Never bare: local branches have drifted from the fork, and assembling
    from them would train the rerere cache on content CI never sees.

    Merges run in manifest order: disjoint features first, the hub last. An
    optional feature that will not merge is dropped and recorded; a required one
    raises. The starting branch is always restored -- on success, on abort, and
    on any unexpected git failure -- with the single deliberate exception of
    `stop_at` below.

    `stop_at` is training mode: stop at that feature and LEAVE its conflicted
    merge in the worktree, so a human resolves it at the exact intermediate state
    a later run will reproduce. Resolving it anywhere else records a pre-image
    that will not match. It is validated against the manifest first: a typo that
    matches nothing would otherwise run an ordinary full build and report
    success, silently doing the opposite of what was asked.
    """
    names = [f.name for f in manifest.features]
    if stop_at is not None and stop_at not in names:
        raise UnknownStopTarget(
            f"--stop-at '{stop_at}' names no feature in the manifest; "
            f"valid names: {', '.join(names) or '(none)'}"
        )

    # Report every missing branch at once. Custody is a property of the whole
    # manifest, and fixing them one failed run at a time is a slow way to learn
    # that a prune took three branches, not one.
    missing = [
        f"{remote}/{f.branch} (feature '{f.name}')"
        for f in manifest.features
        if not git_ok(repo, "rev-parse", "--verify", f"{remote}/{f.branch}^{{commit}}")
    ]
    if missing:
        raise MissingBranch(
            f"missing from the fork: {', '.join(missing)} — "
            "mirror every manifest branch into the fork before assembling"
        )

    # Resolve the base to a commit before anything moves: the release job proves
    # a tag matches the tree it is about to name, and only the exact commit the
    # features were merged onto can settle that. `head` cannot -- it is the tip
    # after every merge, which is never an upstream commit.
    base = stdout(repo, "rev-parse", f"{base_ref}^{{commit}}")
    start_ref = _current_ref(repo)
    git(repo, "checkout", "-q", "-B", out_branch, base_ref)

    merged: list[str] = []
    dropped: list[DroppedFeature] = []
    restore_start_ref = True
    try:
        for feature in manifest.features:
            ok, conflicts = _try_merge(repo, feature, f"{remote}/{feature.branch}")
            if ok:
                merged.append(feature.name)
                continue
            if stop_at == feature.name:
                # Leave the conflict in the worktree on purpose. This is the one
                # path that must NOT restore the starting branch -- the human
                # needs this exact state.
                restore_start_ref = False
                raise TrainingStop(
                    f"stopped at '{feature.name}' for rerere training; "
                    f"resolve {', '.join(conflicts)} here, commit, then reset --hard HEAD~1"
                )
            git(repo, "merge", "--abort", check=False)
            if feature.required:
                raise RequiredFeatureConflict(
                    f"required feature '{feature.name}' conflicts in: {', '.join(conflicts)}"
                )
            dropped.append(DroppedFeature(name=feature.name, conflicting_files=conflicts))
        # Read HEAD inside the `try`, before the `finally` can move it: a restore
        # that ran first would report the caller's own commit as this build's.
        head = stdout(repo, "rev-parse", "HEAD")
    finally:
        # A `finally`, not an `except RequiredFeatureConflict`: a git call that
        # fails for any other reason strands the worktree on `out_branch`, where
        # `minibwa_dist/` does not exist, and the caller is left holding a repo it
        # never asked to be moved.
        if restore_start_ref:
            git(repo, "checkout", "-q", start_ref)

    return AssemblyResult(base=base, head=head, merged=tuple(merged), dropped=tuple(dropped))


_MB_VERSION = re.compile(r'^#define MB_VERSION "([^"]+)"', re.M)


def set_version(repo: Path, version: str) -> str:
    """Rewrite `minibwa.h`'s MB_VERSION to `version`. Returns the previous value.

    The single place that knows the header's shape. The release workflow used to
    carry a second copy as a `sed` expression with looser matching and no
    existence check, so a change to how the define is written had to be found in
    two places and only one of them had tests -- and an operator-supplied version
    containing `&` or `/` was spliced into `sed`'s replacement as syntax. Python's
    `re.sub` takes the replacement as a callable here for the same reason: nothing
    in `version` is ever interpreted.
    """
    header = repo / "minibwa.h"
    text = header.read_text()
    current = _MB_VERSION.search(text)
    if current is None:
        raise RuntimeError("MB_VERSION not found in minibwa.h")
    header.write_text(_MB_VERSION.sub(lambda _: f'#define MB_VERSION "{version}"', text, count=1))
    return current.group(1)


def stamp_dev_version(repo: Path, out_branch: str) -> str:
    """Rewrite MB_VERSION on `out_branch` so every dist build identifies itself.

    Without this, `dist` -- the default branch, what every clone gets -- builds as
    stock upstream and its SAM @PG VN: tag lies about provenance. The release
    workflow re-stamps the same line with the final version.

    Deliberately leaves HEAD on `out_branch` on return -- a caller that needs
    the original branch back is responsible for restoring it itself.
    """
    # Switch FIRST: every read below must come from out_branch. Reading
    # minibwa.h on the calling branch and writing it back after checkout would
    # clobber whatever the assembled branch holds -- and four manifest features
    # (submem-ablation, soft-clip-penalty, meth-cleanups, alt-liftgroup) edit
    # this header. The byte-identity gate would not catch it: reverting toward
    # stock makes the SAM more identical, not less. It also made the embedded
    # sha identify the tooling checkout rather than the build being stamped.
    git(repo, "checkout", "-q", out_branch)

    current = _MB_VERSION.search((repo / "minibwa.h").read_text())
    if current is None:
        raise RuntimeError("MB_VERSION not found in minibwa.h")
    upstream = current.group(1).split("-")[0]
    sha = stdout(repo, "rev-parse", "--short", "HEAD")
    version = f"{upstream}-nh13.dev+{sha}"

    set_version(repo, version)
    git(repo, "add", "--", "minibwa.h")
    git(
        repo,
        "-c",
        "commit.gpgsign=false",
        "commit",
        "-m",
        f"build: identify this build as {version}",
    )
    return version
