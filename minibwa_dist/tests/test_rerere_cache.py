"""Tests for rerere cache persistence on an orphan branch."""

from pathlib import Path

import pytest

from minibwa_dist.rerere_cache import OrphanOverPublishedCache, restore, save
from minibwa_dist.tests.conftest import merge_and_resolve, run


def test_save_returns_false_when_cache_is_empty(repo: Path) -> None:
    assert save(repo) is False


def test_restore_on_missing_branch_returns_zero(repo: Path) -> None:
    assert restore(repo) == 0


def test_round_trip_preserves_entries(repo: Path) -> None:
    """A saved cache must come back byte-identical on a clean clone."""
    rr = repo / ".git" / "rr-cache" / "deadbeef"
    rr.mkdir(parents=True)
    (rr / "preimage").write_text("<<<<<<<\nours\n=======\ntheirs\n>>>>>>>\n")
    (rr / "postimage").write_text("resolved\n")

    assert save(repo) is True
    assert run(repo, "rev-parse", "--verify", "rerere-cache")

    # Simulate a fresh runner: blow the local cache away.
    (rr / "preimage").unlink()
    (rr / "postimage").unlink()
    rr.rmdir()

    assert restore(repo) == 1
    assert (rr / "postimage").read_text() == "resolved\n"


def test_save_is_idempotent(repo: Path) -> None:
    rr = repo / ".git" / "rr-cache" / "cafe"
    rr.mkdir(parents=True)
    (rr / "postimage").write_text("x\n")
    assert save(repo) is True
    assert save(repo) is False, "unchanged cache must not create an empty commit"


def test_save_refuses_to_orphan_over_a_published_cache(repo: Path) -> None:
    """Regression: a swallowed fetch is indistinguishable from a first run.

    When the cache branch failed to fetch, `save()` saw no local branch, built a
    fresh orphan holding only this run's entry, and the caller published it --
    destroying every accumulated resolution in the tip and in the history. The
    remote-tracking ref is the one local signal that separates "never fetched"
    from "nothing to fetch"; with it present, refuse.
    """
    published = repo / ".git" / "rr-cache" / "aaaa1111"
    published.mkdir(parents=True)
    (published / "postimage").write_text("first\n")
    assert save(repo) is True
    run(repo, "update-ref", "refs/remotes/origin/rerere-cache", "rerere-cache")

    # Simulate the swallowed fetch: the branch is published, but this clone
    # never got it into refs/heads.
    run(repo, "update-ref", "-d", "refs/heads/rerere-cache")
    assert restore(repo) == 0, "the precondition is a cache that did not restore"

    fresh = repo / ".git" / "rr-cache" / "cccc3333"
    fresh.mkdir(parents=True)
    (fresh / "postimage").write_text("second\n")

    with pytest.raises(OrphanOverPublishedCache, match="never fetched"):
        save(repo)

    assert not run(repo, "for-each-ref", "--format=%(refname)", "refs/heads/rerere-cache"), (
        "no orphan branch may be left behind for a caller to publish"
    )


def test_save_from_linked_worktree_finds_common_gitdir_cache(repo: Path, tmp_path: Path) -> None:
    """rerere writes its cache under the *common* gitdir even from a linked
    worktree, so `save()` invoked against that worktree must still find it.

    This is exactly this harness's own shape: `ops-distro` is itself a linked
    worktree, where `--absolute-git-dir` would resolve to the private
    `.git/worktrees/<name>` directory rather than the one rerere actually
    writes to.
    """
    run(repo, "checkout", "-b", "feature")
    (repo / "src.c").write_text("line1\nCHANGED-BY-FEATURE\n")
    run(repo, "add", "src.c")
    run(repo, "-c", "commit.gpgsign=false", "commit", "-m", "feature change")

    run(repo, "checkout", "master")
    (repo / "src.c").write_text("line1\nCHANGED-BY-MASTER\n")
    run(repo, "add", "src.c")
    run(repo, "-c", "commit.gpgsign=false", "commit", "-m", "master change")

    worktree = tmp_path / "linked-worktree"
    run(repo, "worktree", "add", "--detach", str(worktree), "master")

    merge_and_resolve(worktree, "feature", "src.c", "line1\nRESOLVED\n")

    assert save(worktree) is True, "save() from a linked worktree found no cache"
    assert run(repo, "rev-parse", "--verify", "rerere-cache")
    assert restore(repo) >= 1
