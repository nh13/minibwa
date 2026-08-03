"""Tests for the merge-assembly engine."""

from pathlib import Path

import pytest

from minibwa_dist.assemble import (
    MergeFailed,
    MissingBranch,
    RequiredFeatureConflict,
    TrainingStop,
    UnknownStopTarget,
    assemble,
    set_version,
    stamp_dev_version,
)
from minibwa_dist.manifest import Feature, Manifest, Upstream
from minibwa_dist.tests.conftest import (
    branch_touching,
    commit_file,
    conflicting_branch,
    merge_and_resolve,
    run,
)


def _feature(name: str, branch: str, *, required: bool = False) -> Feature:
    return Feature(
        name=name,
        branch=branch,
        required=required,
        output="identical",
        summary=name,
        upstream=Upstream(status="unsubmitted"),
    )


def test_merges_disjoint_features(repo: Path) -> None:
    branch_touching(repo, "feat/a", "a.c", "a\n")
    branch_touching(repo, "feat/b", "b.c", "b\n")
    manifest = Manifest(features=(_feature("a", "feat/a"), _feature("b", "feat/b")), withdrawn=())

    result = assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert result.merged == ("a", "b")
    assert result.dropped == ()
    run(repo, "checkout", "-q", "dist-next")
    assert (repo / "a.c").exists() and (repo / "b.c").exists()


def test_resolves_through_the_remote_not_a_stale_local_branch(repo: Path) -> None:
    """The C6 regression: a drifted local branch must never reach the build."""
    branch_touching(repo, "feat/a", "a.c", "fork version\n")
    run(repo, "checkout", "-q", "feat/a")
    commit_file(repo, "a.c", "STALE LOCAL\n", "local drift the fork never saw")
    run(repo, "checkout", "-q", "master")

    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())
    assemble(repo, manifest, "master", "dist-next", remote="fork")

    run(repo, "checkout", "-q", "dist-next")
    assert (repo / "a.c").read_text() == "fork version\n"


def test_restores_the_starting_branch(repo: Path) -> None:
    branch_touching(repo, "feat/a", "a.c", "a\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())
    assemble(repo, manifest, "master", "dist-next", remote="fork")
    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master"


def test_drops_optional_feature_on_conflict(repo: Path) -> None:
    """An optional feature that will not merge is dropped; the build continues."""
    conflicting_branch(repo, "feat/x")
    branch_touching(repo, "feat/ok", "ok.c", "ok\n")

    manifest = Manifest(features=(_feature("x", "feat/x"), _feature("ok", "feat/ok")), withdrawn=())
    result = assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert result.merged == ("ok",)
    assert [d.name for d in result.dropped] == ["x"]
    assert result.dropped[0].conflicting_files == ("src.c",)
    run(repo, "checkout", "-q", "dist-next")
    assert (repo / "ok.c").exists(), "a dropped feature must not abort the build"


def test_required_feature_conflict_aborts_and_restores(repo: Path) -> None:
    conflicting_branch(repo, "feat/req")
    manifest = Manifest(features=(_feature("req", "feat/req", required=True),), withdrawn=())
    with pytest.raises(RequiredFeatureConflict, match="req"):
        assemble(repo, manifest, "master", "dist-next", remote="fork")
    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master", (
        "an abort must not strand the worktree on dist-next, where minibwa_dist/ does not exist"
    )


def test_missing_branch_aborts(repo: Path) -> None:
    """Branch custody: a vanished feature branch is a hard error, never a silent drop."""
    manifest = Manifest(features=(_feature("gone", "feat/gone"),), withdrawn=())
    with pytest.raises(MissingBranch, match="fork/feat/gone"):
        assemble(repo, manifest, "master", "dist-next", remote="fork")


def test_merges_are_unsigned_even_when_gpgsign_is_on(repo: Path) -> None:
    """git merge honours commit.gpgsign; an unavailable key would abort the build."""
    run(repo, "config", "commit.gpgsign", "true")
    run(repo, "config", "gpg.program", "/bin/false")  # any signing attempt fails
    branch_touching(repo, "feat/a", "a.c", "a\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())
    assert assemble(repo, manifest, "master", "dist-next", remote="fork").merged == ("a",)


def test_stop_at_leaves_the_conflict_for_a_human(repo: Path) -> None:
    """Training mode: rerere must learn the pre-image at the exact intermediate state."""
    branch_touching(repo, "feat/first", "first.c", "first\n")
    conflicting_branch(repo, "feat/x")
    manifest = Manifest(
        features=(_feature("first", "feat/first"), _feature("x", "feat/x")), withdrawn=()
    )
    with pytest.raises(RequiredFeatureConflict, match="stopped at 'x'"):
        assemble(repo, manifest, "master", "dist-next", remote="fork", stop_at="x")
    # The conflicted merge is left in place, on top of everything ordered before it.
    assert run(repo, "diff", "--name-only", "--diff-filter=U") == "src.c"
    assert (repo / "first.c").exists()


def test_rerere_replays_a_recorded_resolution(repo: Path) -> None:
    """The load-bearing property of the whole design: a resolved conflict stays resolved.

    If this breaks, every sync costs a human the same conflicts it cost last time.
    """
    run(repo, "config", "rerere.enabled", "true")
    run(repo, "config", "rerere.autoUpdate", "true")
    run(repo, "checkout", "-q", "-b", "feat/x", "master")
    commit_file(repo, "src.c", "line1\nX\n", "x edits line2")
    run(repo, "checkout", "-q", "master")
    commit_file(repo, "src.c", "line1\nUPSTREAM\n", "upstream edits line2")

    run(repo, "update-ref", "refs/remotes/fork/feat/x", "feat/x")
    manifest = Manifest(features=(_feature("x", "feat/x"),), withdrawn=())

    # First pass: the conflict is novel, so the optional feature drops --
    # but rerere records the pre-image on the way past. (Verified against
    # git 2.50.1: the pre-image survives `git merge --abort`.)
    first = assemble(repo, manifest, "master", "dist-next", remote="fork")
    assert [d.name for d in first.dropped] == ["x"]

    # Teach rerere the resolution once, by hand, exactly as a maintainer would.
    run(repo, "checkout", "-q", "master")
    merge_and_resolve(repo, "fork/feat/x", "src.c", "line1\nRESOLVED\n")
    run(repo, "reset", "--hard", "HEAD~1")  # drop the merge; rerere keeps the resolution

    # Second pass: same pre-image, so rerere replays it and the feature merges.
    second = assemble(repo, manifest, "master", "dist-next", remote="fork")
    assert second.merged == ("x",)
    assert second.dropped == ()


def test_rerere_replays_despite_local_zdiff3_conflict_style(repo: Path) -> None:
    """S1 regression: `_try_merge` must not depend on ambient rerere/conflict config.

    Deliberately does NOT set `rerere.enabled`/`rerere.autoUpdate` via `git
    config` anywhere in this test -- only `merge.conflictStyle = zdiff3`, the
    setting a maintainer bootstrapping the cache locally (via `--stop-at`)
    commonly has set globally. `rerere.enabled` defaults to false until
    `.git/rr-cache` already exists (a fresh repo's chicken-and-egg problem), so
    if `_try_merge` relies on ambient config for either rerere knob, training
    records nothing here and the second assembly drops the feature again.

    The `merge.conflictStyle=merge` pin is checked separately and for what it
    actually does. Measured on git 2.50.1, rerere normalises its pre-image, so
    replay survives a train/replay style mismatch and this test would still pass
    with that pin removed -- the docstring used to claim otherwise. What the pin
    genuinely controls is the markers a human is handed during training, so that
    is what is asserted: no `|||||||` base section despite the ambient zdiff3.
    """
    run(repo, "checkout", "-q", "-b", "feat/x", "master")
    commit_file(repo, "src.c", "line1\nX\n", "x edits line2")
    run(repo, "checkout", "-q", "master")
    commit_file(repo, "src.c", "line1\nUPSTREAM\n", "upstream edits line2")
    run(repo, "update-ref", "refs/remotes/fork/feat/x", "feat/x")
    manifest = Manifest(features=(_feature("x", "feat/x"),), withdrawn=())

    # Simulate a maintainer bootstrapping the cache with their own (zdiff3)
    # global preference set locally, via `--stop-at` training mode.
    run(repo, "config", "merge.conflictStyle", "zdiff3")
    with pytest.raises(RequiredFeatureConflict, match="stopped at 'x'"):
        assemble(repo, manifest, "master", "dist-next", remote="fork", stop_at="x")

    conflicted = (repo / "src.c").read_text()
    assert "<<<<<<<" in conflicted, "training must leave the conflict for a human"
    assert "|||||||" not in conflicted, (
        "the pinned `merge` conflict style must win over the maintainer's ambient zdiff3, "
        "so the markers resolved by hand are the ones CI would have produced"
    )

    (repo / "src.c").write_text("line1\nRESOLVED\n")
    run(repo, "add", "src.c")
    run(repo, "-c", "commit.gpgsign=false", "commit", "--no-edit")
    run(repo, "reset", "--hard", "HEAD~1")  # drop the merge; rerere keeps the resolution
    run(repo, "checkout", "-q", "master")

    second = assemble(repo, manifest, "master", "dist-next", remote="fork")
    assert second.merged == ("x",)
    assert second.dropped == ()


def test_the_training_stop_is_its_own_exception_type(repo: Path) -> None:
    """The two RequiredFeatureConflict messages used to be told apart by prefix.

    A subclass keeps `except RequiredFeatureConflict` working while letting the
    CLI dispatch on type, so editing either message cannot silently turn a
    failed build into a reported-successful training stop.
    """
    conflicting_branch(repo, "feat/x")
    manifest = Manifest(features=(_feature("x", "feat/x"),), withdrawn=())
    with pytest.raises(TrainingStop):
        assemble(repo, manifest, "master", "dist-next", remote="fork", stop_at="x")


def test_required_conflict_restores_even_with_stop_at_set(repo: Path) -> None:
    """The restore must key on which failure fired, not on whether --stop-at was passed.

    A maintainer training a later feature whose earlier required feature
    conflicts was left on dist-next -- the exact state
    test_required_feature_conflict_aborts_and_restores calls unacceptable.
    """
    conflicting_branch(repo, "feat/req")
    branch_touching(repo, "feat/later", "later.c", "later\n")
    manifest = Manifest(
        features=(_feature("req", "feat/req", required=True), _feature("later", "feat/later")),
        withdrawn=(),
    )

    with pytest.raises(RequiredFeatureConflict, match="required feature 'req'"):
        assemble(repo, manifest, "master", "dist-next", remote="fork", stop_at="later")

    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master"


def test_an_unexpected_git_failure_still_restores(
    repo: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Anything escaping the merge loop must still hand the worktree back.

    The restore used to live in `except RequiredFeatureConflict`, so any other
    exception -- a git call that failed for a reason nobody enumerated -- left
    the caller parked on dist-next, where `minibwa_dist/` does not exist.
    """
    branch_touching(repo, "feat/a", "a.c", "a\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())

    def boom(*_args: object, **_kwargs: object) -> tuple[bool, tuple[str, ...]]:
        raise RuntimeError("git fell over")

    monkeypatch.setattr("minibwa_dist.assemble._try_merge", boom)
    with pytest.raises(RuntimeError, match="git fell over"):
        assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master"


def test_a_detached_head_is_restored_to_its_commit(repo: Path) -> None:
    """`rev-parse --abbrev-ref HEAD` answers "HEAD" when detached, so restoring
    it was a no-op that silently left the caller on the assembled branch.
    """
    branch_touching(repo, "feat/a", "a.c", "a\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())
    run(repo, "checkout", "-q", "--detach", "master")
    start = run(repo, "rev-parse", "HEAD")

    assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert run(repo, "rev-parse", "HEAD") == start
    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "HEAD", "must still be detached"
    assert not (repo / "a.c").exists(), "the assembled tree must not be left in the worktree"


def test_a_merge_that_never_started_is_not_read_as_a_rerere_hit(repo: Path) -> None:
    """ "Nonzero exit, nothing unmerged" also describes a merge git refused to begin.

    An untracked file the merge would overwrite produces exactly that signature.
    Committing on the strength of it ran `git commit --no-edit` outside a merge,
    which failed with git's own explanation swallowed by the wrapper.
    """
    branch_touching(repo, "feat/a", "a.c", "from the branch\n")
    (repo / "a.c").write_text("untracked, and in the way\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())

    with pytest.raises(MergeFailed, match="never started"):
        assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master"


def test_stop_at_must_name_a_feature(repo: Path) -> None:
    """A typo used to run an ordinary full build and report success."""
    conflicting_branch(repo, "feat/x")
    manifest = Manifest(features=(_feature("x", "feat/x"),), withdrawn=())

    with pytest.raises(UnknownStopTarget, match="names no feature"):
        assemble(repo, manifest, "master", "dist-next", remote="fork", stop_at="x-typo")


def test_every_missing_branch_is_named_at_once(repo: Path) -> None:
    """Custody is a property of the manifest; learning it one run at a time is slow."""
    branch_touching(repo, "feat/here", "here.c", "here\n")
    manifest = Manifest(
        features=(
            _feature("gone", "feat/gone"),
            _feature("here", "feat/here"),
            _feature("also-gone", "feat/also-gone"),
        ),
        withdrawn=(),
    )
    with pytest.raises(MissingBranch) as excinfo:
        assemble(repo, manifest, "master", "dist-next", remote="fork")

    message = str(excinfo.value)
    assert "fork/feat/gone" in message and "fork/feat/also-gone" in message
    assert "feat/here" not in message


def test_the_result_records_the_commit_the_build_was_assembled_from(repo: Path) -> None:
    """`head` is the tip after every merge, so it can never identify the upstream
    commit a release claims to mirror. Only the base can, and the release job
    compares it to the tag it was asked to publish.
    """
    branch_touching(repo, "feat/a", "a.c", "a\n")
    manifest = Manifest(features=(_feature("a", "feat/a"),), withdrawn=())
    expected = run(repo, "rev-parse", "master")

    result = assemble(repo, manifest, "master", "dist-next", remote="fork")

    assert result.base == expected
    assert result.head != expected, "the assembly tip is not the base it was built on"


def test_stamp_dev_version_rewrites_only_the_version_line(repo: Path) -> None:
    """A synthetic header, never the real one -- see the module docstring."""
    header = '#include <foo.h>\n\n#define MB_VERSION "0.6-r416"\n#define OTHER 1\n'
    commit_file(repo, "minibwa.h", header, "add header")
    sha = run(repo, "rev-parse", "--short", "HEAD")

    version = stamp_dev_version(repo, "master")

    assert version == f"0.6-nh13.dev+{sha}"
    assert (repo / "minibwa.h").read_text() == header.replace('"0.6-r416"', f'"{version}"'), (
        "everything but the MB_VERSION line must come through byte-identical"
    )


def test_stamp_dev_version_does_not_compound_on_a_second_run(repo: Path) -> None:
    """The already-stamped case must strip back to the bare upstream version first."""
    header = '#define MB_VERSION "0.6-r416"\n'
    commit_file(repo, "minibwa.h", header, "add header")

    first = stamp_dev_version(repo, "master")
    second = stamp_dev_version(repo, "master")

    assert first.startswith("0.6-nh13.dev+") and first.count("-nh13.dev+") == 1
    assert second.startswith("0.6-nh13.dev+") and second.count("-nh13.dev+") == 1
    assert (repo / "minibwa.h").read_text() == f'#define MB_VERSION "{second}"\n'


def test_stamp_dev_version_raises_when_mb_version_is_absent(repo: Path) -> None:
    commit_file(repo, "minibwa.h", "#define OTHER 1\n", "add header without a version")
    with pytest.raises(RuntimeError, match="MB_VERSION"):
        stamp_dev_version(repo, "master")


def test_stamp_dev_version_reads_and_writes_only_out_branch(repo: Path) -> None:
    """Regression: must not clobber out_branch's header with the calling branch's copy.

    Four manifest features edit minibwa.h (submem-ablation, soft-clip-penalty,
    meth-cleanups, alt-liftgroup). Reading the header on the calling branch and
    writing it back after `checkout out_branch` would silently revert those
    features' edits on every assembled build -- and the byte-identity gate would
    not catch it, since reverting toward stock makes the SAM MORE identical, not
    less. This also pins the embedded sha to `out_branch`'s HEAD, not the caller's.
    """
    commit_file(repo, "minibwa.h", '#define MB_VERSION "0.6-r416"\n', "add header on master")

    run(repo, "checkout", "-q", "-b", "dist-next", "master")
    dist_header = '#define MB_VERSION "0.6-r416"\n#define FEATURE_MARKER 1\n'
    commit_file(repo, "minibwa.h", dist_header, "a feature edits minibwa.h")
    dist_sha = run(repo, "rev-parse", "--short", "HEAD")
    run(repo, "checkout", "-q", "master")

    version = stamp_dev_version(repo, "dist-next")

    assert version == f"0.6-nh13.dev+{dist_sha}", (
        "the embedded sha must identify out_branch's HEAD, not the calling branch's"
    )
    run(repo, "checkout", "-q", "dist-next")
    text = (repo / "minibwa.h").read_text()
    assert "#define FEATURE_MARKER 1" in text, (
        "a feature's edit to minibwa.h must survive stamping, not revert to the "
        "calling branch's stale copy"
    )
    assert text == dist_header.replace('"0.6-r416"', f'"{version}"')


def test_set_version_never_interprets_the_version_string(repo: Path) -> None:
    """The release workflow used to do this rewrite in `sed`, where an
    operator-supplied `&` re-inserted the whole matched line into the header and
    a `/` ended the substitution early. Nothing in the version is syntax here.
    """
    header = '#define MB_VERSION "0.6-r416"\n#define OTHER 1\n'
    (repo / "minibwa.h").write_text(header)

    previous = set_version(repo, "0.6-nh13.a&b/c\\d")

    assert previous == "0.6-r416"
    assert (repo / "minibwa.h").read_text() == (
        '#define MB_VERSION "0.6-nh13.a&b/c\\d"\n#define OTHER 1\n'
    )


def test_set_version_refuses_a_header_without_the_define(repo: Path) -> None:
    (repo / "minibwa.h").write_text("#define OTHER 1\n")
    with pytest.raises(RuntimeError, match="MB_VERSION"):
        set_version(repo, "0.6-nh13.1")
