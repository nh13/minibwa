"""Tests for the byte-identity output gates."""

from pathlib import Path

import pytest

from minibwa_dist.gates import run_gates, sam_digest
from minibwa_dist.manifest import Feature, Manifest, Upstream


def _sam(path: Path, *lines: str) -> Path:
    path.write_text("".join(line + "\n" for line in lines))
    return path


def test_digest_ignores_pg_and_co(tmp_path: Path) -> None:
    """@PG carries MB_VERSION and argv, which differ by construction."""
    a = _sam(
        tmp_path / "a.sam",
        "@HD\tVN:1.6",
        "@PG\tID:minibwa\tVN:0.6-r416",
        "@CO\tanything",
        "r1\t0\tchrM\t1\t60\t5M",
    )
    b = _sam(
        tmp_path / "b.sam",
        "@HD\tVN:1.6",
        "@PG\tID:minibwa\tVN:0.6-nh13.1",
        "r1\t0\tchrM\t1\t60\t5M",
    )
    assert sam_digest(a) == sam_digest(b)


def test_digest_detects_an_alignment_change(tmp_path: Path) -> None:
    a = _sam(tmp_path / "a.sam", "@HD\tVN:1.6", "r1\t0\tchrM\t1\t60\t5M")
    b = _sam(tmp_path / "b.sam", "@HD\tVN:1.6", "r1\t0\tchrM\t1\t42\t5M")
    assert sam_digest(a) != sam_digest(b), "a MAPQ change must fail the gate"


def test_digest_detects_a_header_change(tmp_path: Path) -> None:
    """@SQ/@RG/@HD are real output; only @PG and @CO are exempt."""
    a = _sam(tmp_path / "a.sam", "@SQ\tSN:chrM\tLN:16569", "r1\t0\tchrM\t1\t60\t5M")
    b = _sam(tmp_path / "b.sam", "@SQ\tSN:chrM\tLN:16570", "r1\t0\tchrM\t1\t60\t5M")
    assert sam_digest(a) != sam_digest(b)


def _stub_aligner(path: Path, sam: str) -> Path:
    """A `minibwa` stand-in: `index` succeeds, `map` prints a fixed SAM.

    The comparison step -- the single property this module exists for -- was
    unreachable from the tests, because reaching it needs two aligners. It does
    not need two *real* ones.
    """
    path.write_text(
        '#!/bin/sh\nif [ "$1" = "index" ]; then exit 0; fi\ncat <<\'SAM\'\n' + sam + "SAM\n"
    )
    path.chmod(0o755)
    return path


def _fixtures(tmp_path: Path) -> Path:
    directory = tmp_path / "fixtures"
    directory.mkdir()
    for name in ("chrM-human.fa.gz", "chrM-read_1.fa.gz", "chrM-read_2.fa.gz"):
        (directory / name).write_bytes(b"")
    return directory


def _identical_feature(name: str) -> Feature:
    return Feature(
        name=name,
        branch=f"feat/{name}",
        required=False,
        output="identical",
        summary=name,
        upstream=Upstream(status="unsubmitted"),
    )


def test_the_gate_passes_identical_output_and_fails_a_difference(tmp_path: Path) -> None:
    """Hard-wiring `passed = True` used to leave the whole suite green: nothing
    drove `run_gates` far enough to compare two outputs at all.
    """
    fixtures = _fixtures(tmp_path)
    sam = "@HD\tVN:1.6\nr1\t0\tchrM\t1\t60\t5M\n"
    stock = _stub_aligner(tmp_path / "stock", sam)
    manifest = Manifest(features=(_identical_feature("a"),), withdrawn=())

    same = _stub_aligner(tmp_path / "same", sam)
    [result] = run_gates(same, stock, fixtures, manifest, tmp_path / "pass")
    assert result.passed is True
    assert result.detail == "1 feature(s) covered: a"

    differs = _stub_aligner(tmp_path / "differs", sam.replace("\t60\t", "\t42\t"))
    [result] = run_gates(differs, stock, fixtures, manifest, tmp_path / "fail")
    assert result.passed is False
    assert "SAM differs" in result.detail


def test_a_relative_binary_path_still_resolves(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Both workflows pass `Path("./minibwa")`, which pathlib normalises to
    `minibwa` -- a slash-free string subprocess looks up on PATH, not in cwd.
    """
    fixtures = _fixtures(tmp_path)
    _stub_aligner(tmp_path / "minibwa", "@HD\tVN:1.6\n")
    monkeypatch.chdir(tmp_path)

    [result] = run_gates(
        Path("./minibwa"),
        Path("./minibwa"),
        fixtures,
        Manifest(features=(), withdrawn=()),
        tmp_path / "work",
    )
    assert result.passed is True, result.detail


def test_coverage_names_only_the_features_in_this_build(tmp_path: Path) -> None:
    """An optional feature drops out routinely; claiming the gate covered it
    credits code the binary does not contain.
    """
    fixtures = _fixtures(tmp_path)
    binary = _stub_aligner(tmp_path / "minibwa", "@HD\tVN:1.6\n")
    manifest = Manifest(
        features=(_identical_feature("kept"), _identical_feature("dropped")), withdrawn=()
    )

    [result] = run_gates(binary, binary, fixtures, manifest, tmp_path / "work", merged=("kept",))

    assert result.passed is True
    assert result.detail == "1 feature(s) covered: kept; 1 not in this build: dropped"


def test_a_dropped_changes_output_feature_does_not_block_the_build(tmp_path: Path) -> None:
    """The always-fail gate exists to stop unmeasured output changes shipping. A
    feature that was dropped is not shipping, so it has nothing to block.
    """
    manifest = Manifest(
        features=(
            Feature(
                name="risky",
                branch="feat/risky",
                required=False,
                output="changes-output",
                summary="x",
                upstream=Upstream(status="unsubmitted"),
            ),
        ),
        withdrawn=(),
    )
    results = run_gates(
        Path("/nonexistent"), Path("/nonexistent"), tmp_path, manifest, tmp_path, merged=()
    )
    assert not [r for r in results if r.name == "changes-output:risky"]


def test_changes_output_features_always_fail(tmp_path: Path) -> None:
    """Design section 8: a changes-output feature is blocked until measured."""
    manifest = Manifest(
        features=(
            Feature(
                name="risky",
                branch="feat/risky",
                required=False,
                output="changes-output",
                summary="x",
                upstream=Upstream(status="unsubmitted"),
            ),
        ),
        withdrawn=(),
    )
    results = run_gates(Path("/nonexistent"), Path("/nonexistent"), tmp_path, manifest, tmp_path)
    blocked = [r for r in results if r.name == "changes-output:risky"]
    assert blocked and blocked[0].passed is False


def test_binary_failure_is_a_gate_failure_not_a_traceback(tmp_path: Path) -> None:
    """A crashing aligner must report FAIL, not blow up the workflow."""
    results = run_gates(
        Path("/nonexistent"),
        Path("/nonexistent"),
        tmp_path,
        Manifest(features=(), withdrawn=()),
        tmp_path,
    )
    assert any(not r.passed and "could not run" in r.detail for r in results)
