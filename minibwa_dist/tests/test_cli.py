"""Tests for the CLI seam (`minibwa_dist.cli.main`), using synthetic repos.

No `gh` calls and no network here -- these exercise argument parsing and the
glue between `assemble()`/`stamp_dev_version()` and the process exit code, the
same way `distro-sync.yml` invokes `python -m minibwa_dist.cli sync`.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import minibwa_dist.cli as cli_module
from minibwa_dist.cli import main
from minibwa_dist.manifest import load_manifest
from minibwa_dist.render import BLOCK_BEGIN
from minibwa_dist.tests.conftest import (
    branch_touching,
    commit_file,
    conflicting_branch,
    run,
    stub_aligner,
)

_ONE_FEATURE = """
[[feature]]
name = "{name}"
branch = "feat/{name}"
required = {required}
output = "identical"
summary = "{name}"
upstream = {{ status = "unsubmitted" }}
"""


def _write_manifest(tmp_path: Path, features_toml: str) -> Path:
    path = tmp_path / "features.toml"
    path.write_text(features_toml)
    return path


def _feature_toml(name: str, *, required: bool = False) -> str:
    return _ONE_FEATURE.format(name=name, required="true" if required else "false")


def _sync_argv(manifest: Path, repo: Path, *extra: str) -> list[str]:
    return [
        "--manifest",
        str(manifest),
        "sync",
        "--repo",
        str(repo),
        "--base",
        "master",
        "--out",
        "dist-next",
        "--remote",
        "fork",
        *extra,
    ]


def test_sync_cold_assembly_merges_disjoint_features(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A fresh repo, no rerere cache: the ordinary cold-start path."""
    branch_touching(repo, "feat/a", "a.c", "a\n")
    branch_touching(repo, "feat/b", "b.c", "b\n")
    base_sha = run(repo, "rev-parse", "master")
    manifest = _write_manifest(
        tmp_path,
        """
[[feature]]
name = "a"
branch = "feat/a"
required = false
output = "identical"
summary = "a"
upstream = { status = "unsubmitted" }

[[feature]]
name = "b"
branch = "feat/b"
required = false
output = "identical"
summary = "b"
upstream = { status = "unsubmitted" }
""",
    )

    rc = main(_sync_argv(manifest, repo))

    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["merged"] == ["a", "b"]
    assert out["dropped"] == []
    assert out["version"] is None
    assert out["base"] == base_sha, (
        "the release job proves a tag matches the tree it names by reading this field "
        "back out of dist-manifest.json"
    )
    run(repo, "checkout", "-q", "dist-next")
    assert (repo / "a.c").exists() and (repo / "b.c").exists()


def test_sync_drops_an_optional_feature_and_reports_it(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(
        tmp_path,
        """
[[feature]]
name = "x"
branch = "feat/x"
required = false
output = "identical"
summary = "x"
upstream = { status = "unsubmitted" }
""",
    )

    rc = main(_sync_argv(manifest, repo))

    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["merged"] == []
    assert out["dropped"] == [{"name": "x", "files": ["src.c"]}]


def test_sync_with_stamp_produces_the_expected_version_string(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    commit_file(repo, "minibwa.h", '#define MB_VERSION "0.6-r416"\n', "add header")
    expected_sha = run(repo, "rev-parse", "--short", "master")
    manifest = _write_manifest(tmp_path, "")

    rc = main(_sync_argv(manifest, repo, "--stamp"))

    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["version"] == f"0.6-nh13.dev+{expected_sha}"
    run(repo, "checkout", "-q", "dist-next")
    assert (repo / "minibwa.h").read_text() == f'#define MB_VERSION "{out["version"]}"\n'


def test_sync_stop_at_exits_zero_and_explains_training_mode(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """S10: the expected `--stop-at` outcome must not read as a build failure."""
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(
        tmp_path,
        """
[[feature]]
name = "x"
branch = "feat/x"
required = false
output = "identical"
summary = "x"
upstream = { status = "unsubmitted" }
""",
    )

    rc = main(_sync_argv(manifest, repo, "--stop-at", "x"))

    assert rc == 0
    err = capsys.readouterr().err
    assert "stopped at 'x' for rerere training" in err
    assert "ASSEMBLY FAILED" not in err


def test_sync_required_conflict_still_fails_even_with_stop_at_set(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A genuine required-feature conflict must still exit 1, even when
    `--stop-at` names a different (later, never-reached) feature.

    `later` is a real manifest entry, ordered after the required feature that
    fails: `--stop-at` is now validated against the manifest up front, so a name
    that matches nothing tests the typo path (below), not this one.
    """
    conflicting_branch(repo, "feat/req")
    branch_touching(repo, "feat/later", "later.c", "later\n")
    manifest = _write_manifest(
        tmp_path, _feature_toml("req", required=True) + _feature_toml("later")
    )

    rc = main(_sync_argv(manifest, repo, "--stop-at", "later"))

    assert rc == 1
    err = capsys.readouterr().err
    assert "ASSEMBLY FAILED" in err
    assert "required feature 'req'" in err
    assert run(repo, "rev-parse", "--abbrev-ref", "HEAD") == "master", (
        "a failed build must hand the worktree back, whether or not --stop-at was set"
    )


def test_sync_rejects_a_stop_at_that_names_no_feature(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A typo used to run an ordinary full build, drop the feature the maintainer
    meant to train on, and exit 0 -- the opposite of what was asked, silently.
    """
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(tmp_path, _feature_toml("x"))

    rc = main(_sync_argv(manifest, repo, "--stop-at", "x-typo"))

    assert rc == 1
    err = capsys.readouterr().err
    assert "ASSEMBLY FAILED" in err
    assert "names no feature" in err
    assert "valid names: x" in err


def test_sync_missing_branch_exits_one(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Branch custody at the CLI seam: distro-sync.yml relies on this exit code
    to stop before it builds, gates and pushes a tree assembled without it.
    """
    manifest = _write_manifest(tmp_path, _feature_toml("gone"))

    rc = main(_sync_argv(manifest, repo))

    assert rc == 1
    assert "ASSEMBLY FAILED" in capsys.readouterr().err


def test_sync_training_mode_names_the_remote_to_publish_the_cache_to(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """`save()` runs at the moment of the conflict, so it captures a pre-image and
    nothing else: the resolution only reaches CI if the maintainer re-runs sync
    and pushes. And `origin` is UPSTREAM in a local checkout, so the follow-up
    has to name the remote the caller actually chose.
    """
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(tmp_path, _feature_toml("x"))

    assert main(_sync_argv(manifest, repo, "--stop-at", "x")) == 0

    err = capsys.readouterr().err
    assert "git push fork rerere-cache" in err
    assert "Re-run this same sync" in err


def test_sync_restores_and_then_replays_the_rerere_cache_branch(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The cache wiring in `_sync`: save-on-the-way-out, restore-on-the-way-in.

    Deleting either call is invisible to a suite that only tests
    `rerere_cache.save`/`restore` in isolation -- the nightly sync still exits 0
    and still emits JSON, it just re-hits every conflict and drops the feature.
    So drive it the way a maintainer does: train, resolve, sync once to persist
    the completed resolution onto the cache branch, wipe the local `rr-cache`
    the way a fresh runner has it, and sync again.
    """
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(tmp_path, _feature_toml("x"))

    assert main(_sync_argv(manifest, repo, "--stop-at", "x")) == 0
    (repo / "src.c").write_text("line1\nRESOLVED\n")
    run(repo, "add", "src.c")
    run(repo, "-c", "commit.gpgsign=false", "commit", "--no-edit")
    run(repo, "reset", "--hard", "HEAD~1")
    run(repo, "checkout", "-q", "master")

    # A normal sync: this is the run whose `save()` sees a complete resolution.
    assert main(_sync_argv(manifest, repo)) == 0
    assert json.loads(capsys.readouterr().out)["merged"] == ["x"], "rerere should replay locally"
    assert run(repo, "rev-parse", "--verify", "refs/heads/rerere-cache"), "save() must have run"

    # Now be a cold runner: the cache branch exists, the local cache does not.
    cache = Path(run(repo, "rev-parse", "--path-format=absolute", "--git-common-dir")) / "rr-cache"
    shutil.rmtree(cache)

    assert main(_sync_argv(manifest, repo)) == 0
    out = capsys.readouterr()
    assert json.loads(out.out)["merged"] == ["x"], (
        "with no local rr-cache the feature can only merge if restore() repopulated it"
    )
    assert "rerere: restored 1 resolution(s)" in out.err


def test_sync_reports_orphan_over_published_cache_instead_of_a_traceback(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A fresh `git clone` is exactly the state `OrphanOverPublishedCache`
    exists to catch: the cache branch is published upstream (a remote-tracking
    ref exists), but this clone's `git fetch` never landed it locally as
    `refs/heads/rerere-cache`. `_sync` calls `rerere_cache.save()`
    unconditionally from its `finally`, so that exception used to propagate
    straight out of `main()` as an unhandled traceback -- discarding the JSON
    of an assembly that had actually succeeded. `main()` must turn it into an
    actionable message and exit 1 instead.
    """
    rr = repo / ".git" / "rr-cache" / "aaaa1111"
    rr.mkdir(parents=True)
    (rr / "postimage").write_text("some resolution\n")
    tip = run(repo, "rev-parse", "master")
    run(repo, "update-ref", "refs/remotes/fork/rerere-cache", tip)
    manifest = _write_manifest(tmp_path, "")

    rc = main(_sync_argv(manifest, repo))

    assert rc == 1
    err = capsys.readouterr().err
    assert "rerere-cache" in err
    assert "fetch" in err.lower(), "the message must say to fetch the branch, not just fail"


def test_reconcile_write_emits_a_manifest_that_loads(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """What the bot commits, the engine must be able to read.

    `distro-reconcile.yml` runs exactly this, then commits, force-pushes and
    opens a PR -- with no step that tries to load the result. A graduated
    feature written back as `status = "merged"` therefore shipped a manifest
    nothing could parse, and merging that PR wedged sync, reconcile and test.
    """
    manifest = _write_manifest(
        tmp_path,
        """
# --- a heading that introduces the block below ---
[[feature]]
name = "graduate"
branch = "feat/graduate"
required = false
output = "identical"
summary = "upstream took this one"
upstream = { pr = 33, status = "open" }

# --- a heading that must survive the deletion above ---
[[feature]]
name = "stays"
branch = "feat/stays"
required = false
output = "identical"
summary = "still ours"
upstream = { pr = 34, status = "open" }
""",
    )
    states = tmp_path / "prs.json"
    states.write_text(json.dumps({"33": "merged", "34": "closed"}))

    rc = main(["--manifest", str(manifest), "reconcile", "--pr-states", str(states), "--write"])

    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["graduated"] == ["graduate"]

    reloaded = load_manifest(manifest)
    assert [f.name for f in reloaded.features] == ["stays"]
    assert reloaded.features[0].upstream.status == "rejected"
    assert "a heading that must survive the deletion above" in manifest.read_text()


def test_reconcile_write_refuses_a_rewrite_it_cannot_load_back(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """`apply_changes` is tested elsewhere to always hand back loadable text,
    but `_reconcile` must not simply trust that invariant at the one point it
    is about to become a committed file: `distro-reconcile.yml` commits,
    force-pushes and opens a PR with no step that loads the result, so this
    check is the only thing standing between a bad rewrite and a manifest
    nothing can parse landing in git history.

    Simulated by monkeypatching `apply_changes` to hand back malformed text --
    not because today's `apply_changes` is known to produce it, but because
    the CLI's own refusal must not depend on that ever remaining true.
    """
    monkeypatch.setattr(cli_module, "apply_changes", lambda text, changes: "not [valid toml")
    manifest = _write_manifest(
        tmp_path,
        _feature_toml("a").replace(
            'upstream = { status = "unsubmitted" }', 'upstream = { pr = 33, status = "open" }'
        ),
    )
    original = manifest.read_text()
    states = tmp_path / "prs.json"
    states.write_text(json.dumps({"33": "closed"}))

    rc = main(["--manifest", str(manifest), "reconcile", "--pr-states", str(states), "--write"])

    assert rc == 1
    assert "ERROR" in capsys.readouterr().err
    assert manifest.read_text() == original, (
        "a rewrite that cannot be read back must not be written"
    )


def test_reconcile_without_write_leaves_the_manifest_untouched(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The read-only preview a maintainer runs first must stay read-only."""
    text = _feature_toml("a").replace(
        'upstream = { status = "unsubmitted" }', 'upstream = { pr = 33, status = "open" }'
    )
    manifest = _write_manifest(tmp_path, text)
    states = tmp_path / "prs.json"
    states.write_text(json.dumps({"33": "closed"}))

    rc = main(["--manifest", str(manifest), "reconcile", "--pr-states", str(states)])

    assert rc == 0
    assert manifest.read_text() == text, "no --write means no write"
    out = json.loads(capsys.readouterr().out)
    assert out["changes"] == [{"feature": "a", "old": "open", "new": "rejected"}]


def test_notes_render_the_dropped_section(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    manifest = _write_manifest(tmp_path, _feature_toml("a"))
    upstream = tmp_path / "upstream.md"
    upstream.write_text("upstream said things")

    rc = main(
        [
            "--manifest",
            str(manifest),
            "notes",
            "--upstream-notes",
            str(upstream),
            "--upstream-tag",
            "v0.6",
            "--dropped",
            "a",
        ]
    )

    assert rc == 0
    out = capsys.readouterr().out
    assert "upstream said things" in out
    assert "not included" in out.lower() and "`a`" in out


def test_docs_writes_the_readme_block_and_the_graveyard(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, _feature_toml("a"))
    readme = tmp_path / "README.md"
    readme.write_text("# minibwa\n\nUpstream prose.\n")
    graveyard = tmp_path / "GRAVEYARD.md"

    rc = main(
        [
            "--manifest",
            str(manifest),
            "docs",
            "--readme",
            str(readme),
            "--graveyard",
            str(graveyard),
        ]
    )

    assert rc == 0
    assert BLOCK_BEGIN in readme.read_text() and "Upstream prose." in readme.read_text()
    assert graveyard.read_text().startswith("# Graveyard")


def test_stamp_sets_an_exact_version(repo: Path, tmp_path: Path) -> None:
    """The release workflow's version rewrite, no longer a second copy in sed."""
    (repo / "minibwa.h").write_text('#define MB_VERSION "0.6-dev+abc123"\n')
    manifest = _write_manifest(tmp_path, "")

    rc = main(
        ["--manifest", str(manifest), "stamp", "--repo", str(repo), "--version", "0.6-nh13.1"]
    )

    assert rc == 0
    assert (repo / "minibwa.h").read_text() == '#define MB_VERSION "0.6-nh13.1"\n'


_stub_aligner = stub_aligner


def test_gates_subcommand_exits_nonzero_when_the_sam_differs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The gate's verdict is the whole point of the step, and it used to live as
    two hand-copied heredocs no linter or test could reach.
    """
    fixtures = tmp_path / "fixtures"
    fixtures.mkdir()
    for name in ("chrM-human.fa.gz", "chrM-read_1.fa.gz", "chrM-read_2.fa.gz"):
        (fixtures / name).write_bytes(b"")
    stock = _stub_aligner(tmp_path / "stock", "@HD\tVN:1.6\nr1\t0\tchrM\t1\t60\t5M\n")
    same = _stub_aligner(tmp_path / "same", "@HD\tVN:1.6\nr1\t0\tchrM\t1\t60\t5M\n")
    different = _stub_aligner(tmp_path / "different", "@HD\tVN:1.6\nr1\t0\tchrM\t1\t42\t5M\n")
    manifest = _write_manifest(tmp_path, _feature_toml("a"))

    def gates(candidate: Path, workdir: str) -> int:
        return main(
            [
                "--manifest",
                str(manifest),
                "gates",
                "--candidate",
                str(candidate),
                "--stock",
                str(stock),
                "--fixtures",
                str(fixtures),
                "--workdir",
                str(tmp_path / workdir),
                "--repo",
                str(tmp_path),
            ]
        )

    assert gates(same, "pass") == 0
    assert "PASS  default-flags-byte-identity" in capsys.readouterr().out

    assert gates(different, "fail") == 1
    assert "FAIL  default-flags-byte-identity" in capsys.readouterr().out


def test_gates_scope_coverage_to_the_features_this_build_contains(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A dropped feature is not in the binary, so naming it as covered overstates
    what the one trustworthy artifact proved -- and contradicts the drop list the
    same job publishes.
    """
    fixtures = tmp_path / "fixtures"
    fixtures.mkdir()
    for name in ("chrM-human.fa.gz", "chrM-read_1.fa.gz", "chrM-read_2.fa.gz"):
        (fixtures / name).write_bytes(b"")
    binary = _stub_aligner(tmp_path / "aligner", "@HD\tVN:1.6\n")
    manifest = _write_manifest(tmp_path, _feature_toml("kept") + _feature_toml("dropped"))
    assembly = tmp_path / "assembly.json"
    assembly.write_text(json.dumps({"merged": ["kept"], "dropped": [{"name": "dropped"}]}))

    rc = main(
        [
            "--manifest",
            str(manifest),
            "gates",
            "--candidate",
            str(binary),
            "--stock",
            str(binary),
            "--fixtures",
            str(fixtures),
            "--workdir",
            str(tmp_path / "work"),
            "--assembly",
            str(assembly),
            "--repo",
            str(tmp_path),
        ]
    )

    assert rc == 0
    detail = capsys.readouterr().out
    assert "1 feature(s) covered: kept" in detail
    assert "1 not in this build: dropped" in detail


# --- a feature that merged in the previous build and does not now ---


def _record_previous_build(repo: Path, ref: str, merged: list[str]) -> None:
    """Commit a prior build's assembly JSON on `ref`, the way sync itself does."""
    run(repo, "checkout", "-q", "-B", ref, "master")
    commit_file(repo, "dist-manifest.json", json.dumps({"merged": merged, "dropped": []}), "record")
    run(repo, "checkout", "-q", "master")


_ONE_CONFLICTING_FEATURE = """
[[feature]]
name = "x"
branch = "feat/x"
required = false
output = "identical"
summary = "x"
upstream = { status = "unsubmitted" }
"""


def test_sync_fails_when_a_feature_stops_merging(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The drop is a regression: the previous build had it, this one does not.

    Exiting zero here is what let a real drop pass unnoticed -- the build then
    matched `dist`, the caller short-circuited on "nothing to ship", and the
    open sync PR kept advertising the feature.
    """
    conflicting_branch(repo, "feat/x")
    _record_previous_build(repo, "dist-prev", ["x"])
    manifest = _write_manifest(tmp_path, _ONE_CONFLICTING_FEATURE)

    rc = main(_sync_argv(manifest, repo, "--regression-baseline", "dist-prev"))

    assert rc == 1
    captured = capsys.readouterr()
    assert json.loads(captured.out)["regressed"] == ["x"]
    assert "x" in captured.err


def test_sync_still_emits_the_assembly_json_when_it_fails(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The caller redirects stdout to a file and reads it after the step fails."""
    conflicting_branch(repo, "feat/x")
    _record_previous_build(repo, "dist-prev", ["x"])
    manifest = _write_manifest(tmp_path, _ONE_CONFLICTING_FEATURE)

    main(_sync_argv(manifest, repo, "--regression-baseline", "dist-prev"))

    out = json.loads(capsys.readouterr().out)
    assert [d["name"] for d in out["dropped"]] == ["x"]
    assert out["merged"] == []


def test_sync_succeeds_when_a_never_merged_feature_is_dropped_again(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A feature that did not merge last time either is not a regression."""
    conflicting_branch(repo, "feat/x")
    _record_previous_build(repo, "dist-prev", [])
    manifest = _write_manifest(tmp_path, _ONE_CONFLICTING_FEATURE)

    rc = main(_sync_argv(manifest, repo, "--regression-baseline", "dist-prev"))

    assert rc == 0
    assert json.loads(capsys.readouterr().out)["regressed"] == []


def test_sync_succeeds_with_no_baseline_ref(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """First run: nothing to compare against, so nothing can have regressed."""
    conflicting_branch(repo, "feat/x")
    manifest = _write_manifest(tmp_path, _ONE_CONFLICTING_FEATURE)

    rc = main(_sync_argv(manifest, repo, "--regression-baseline", "no-such-branch"))

    assert rc == 0
    assert json.loads(capsys.readouterr().out)["regressed"] == []


def test_sync_without_the_flag_reports_no_regressions(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The flag is opt-in; omitting it must not change the exit code."""
    conflicting_branch(repo, "feat/x")
    _record_previous_build(repo, "dist-prev", ["x"])
    manifest = _write_manifest(tmp_path, _ONE_CONFLICTING_FEATURE)

    rc = main(_sync_argv(manifest, repo))

    assert rc == 0
    assert json.loads(capsys.readouterr().out)["regressed"] == []
