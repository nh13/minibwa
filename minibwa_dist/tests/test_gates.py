"""Tests for the byte-identity output gates."""

from pathlib import Path

import pytest

from minibwa_dist.gates import GateResult, run_feature_suites, run_gates, sam_digest
from minibwa_dist.manifest import Feature, Manifest, Upstream
from minibwa_dist.tests.conftest import stub_aligner


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


_stub_aligner = stub_aligner


def _stub_aligner_varying(path: Path, sam: str, flag: str, other: str) -> Path:
    """A `minibwa` stand-in whose output depends on whether `flag` was passed.

    Models the exact failure the extra `_MODES` entries exist to catch: a build
    that matches stock under default flags and diverges under one upstream flag,
    which a paired default-flags-only comparison passes.

    Emits the same flag line as `_stub_aligner` so the two agree in every mode
    except the one carrying `flag`.
    """
    path.write_text(
        '#!/bin/sh\nif [ "$1" = "index" ]; then exit 0; fi\n'
        f'case " $* " in *" {flag} "*) cat <<\'ALT\'\n{other}ALT\n'
        f";; *) cat <<'SAM'\n{sam}SAM\n;; esac\n"
        'flags=""; n=0\n'
        'for a in "$@"; do case "$a" in -*) flags="$flags $a" ;; *) n=$((n+1)) ;; esac; done\n'
        'printf \'mode\\t%s\\t%s\\n\' "$flags" "$n"\n'
    )
    path.chmod(0o755)
    return path


def _default_gate(results: list[GateResult]) -> GateResult:
    """The `default-flags-byte-identity` result, whichever order gates ran in."""
    [result] = [r for r in results if r.name == "default-flags-byte-identity"]
    return result


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
    results = run_gates(same, stock, fixtures, manifest, tmp_path / "pass")
    assert all(r.passed for r in results), [r.detail for r in results if not r.passed]
    assert [r.name for r in results] == [
        "modes-are-distinct",
        "default-flags-byte-identity",
        "byte-identity:no-unmapped",
        "byte-identity:base-tag",
        "byte-identity:eqx-cigar",
        "byte-identity:single-end",
    ]
    assert _default_gate(results).detail == "1 feature(s) covered: a"

    differs = _stub_aligner(tmp_path / "differs", sam.replace("\t60\t", "\t42\t"))
    results = run_gates(differs, stock, fixtures, manifest, tmp_path / "fail")
    # Every byte-identity mode fails. `modes-are-distinct` still passes: it reads
    # stock alone, so a candidate regression neither triggers nor masks it.
    assert not any(r.passed for r in results if r.name != "modes-are-distinct")
    assert "SAM differs" in _default_gate(results).detail


def test_a_relative_binary_path_still_resolves(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Both workflows pass `Path("./minibwa")`, which pathlib normalises to
    `minibwa` -- a slash-free string subprocess looks up on PATH, not in cwd.
    """
    fixtures = _fixtures(tmp_path)
    _stub_aligner(tmp_path / "minibwa", "@HD\tVN:1.6\n")
    monkeypatch.chdir(tmp_path)

    results = run_gates(
        Path("./minibwa"),
        Path("./minibwa"),
        fixtures,
        Manifest(features=(), withdrawn=()),
        tmp_path / "work",
    )
    assert all(r.passed for r in results), [r.detail for r in results if not r.passed]


def test_coverage_names_only_the_features_in_this_build(tmp_path: Path) -> None:
    """An optional feature drops out routinely; claiming the gate covered it
    credits code the binary does not contain.
    """
    fixtures = _fixtures(tmp_path)
    binary = _stub_aligner(tmp_path / "minibwa", "@HD\tVN:1.6\n")
    manifest = Manifest(
        features=(_identical_feature("kept"), _identical_feature("dropped")), withdrawn=()
    )

    results = run_gates(binary, binary, fixtures, manifest, tmp_path / "work", merged=("kept",))

    result = _default_gate(results)
    assert result.passed is True
    assert result.detail == "1 feature(s) covered: kept; 1 not in this build: dropped"


def test_a_divergence_under_only_one_flag_fails_only_that_mode(tmp_path: Path) -> None:
    """The reason `_MODES` has more than one entry: a build can match stock under
    default flags and diverge under `--eqx` alone. The gate must catch it and name
    which invocation exposed it.
    """
    fixtures = _fixtures(tmp_path)
    sam = "@HD\tVN:1.6\nr1\t0\tchrM\t1\t60\t5M\n"
    stock = _stub_aligner(tmp_path / "stock", sam)
    candidate = _stub_aligner_varying(tmp_path / "cand", sam, "--eqx", sam.replace("5M", "3=1X1="))

    results = run_gates(
        candidate,
        stock,
        fixtures,
        Manifest(features=(_identical_feature("a"),), withdrawn=()),
        tmp_path / "work",
    )

    failed = [r for r in results if not r.passed]
    assert [r.name for r in failed] == ["byte-identity:eqx-cigar"]
    assert "--eqx" in failed[0].detail
    assert _default_gate(results).passed is True, "the old gate would have shipped this"


def test_a_redundant_mode_fails_the_distinctness_gate(tmp_path: Path) -> None:
    """`-a` was nearly added as a mode despite being a no-op in `map`. A mode that
    does not change stock's output inflates the coverage claim without adding
    evidence, so it must fail rather than quietly pass five times over.

    Uses a flag-INSENSITIVE stub -- every mode gets byte-identical output, which
    is exactly what a list of no-op flags would produce against a real aligner.
    """
    fixtures = _fixtures(tmp_path)
    flat = tmp_path / "flat"
    flat.write_text(
        '#!/bin/sh\nif [ "$1" = "index" ]; then exit 0; fi\necho "r1\t0\tchrM\t1\t60\t5M"\n'
    )
    flat.chmod(0o755)

    results = run_gates(
        flat, flat, fixtures, Manifest(features=(), withdrawn=()), tmp_path / "work"
    )

    [distinct] = [r for r in results if r.name == "modes-are-distinct"]
    assert distinct.passed is False
    assert "redundant mode(s)" in distinct.detail
    # The byte-identity gates still pass -- the binary matches itself. Only the
    # claim that five modes were exercised is false, which is the point.
    assert all(r.passed for r in results if r.name != "modes-are-distinct")


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


# --- feature test suites: the positive-path coverage byte-identity cannot give ---


def _suite(directory: Path, name: str, *, exit_code: int, message: str = "") -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(f"#!/bin/sh\necho '{message}'\nexit {exit_code}\n")
    path.chmod(0o755)
    return path


def _with_tests(name: str, tests: str | None) -> Feature:
    return Feature(
        name=name,
        branch=f"feat/{name}",
        required=False,
        output="identical",
        summary=name,
        upstream=Upstream(status="unsubmitted"),
        tests=tests,
    )


def test_a_merged_features_suites_run_and_pass(tmp_path: Path) -> None:
    _suite(tmp_path / "test/x", "test-a.sh", exit_code=0)
    _suite(tmp_path / "test/x", "test-b.sh", exit_code=0)
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())

    results = run_feature_suites(tmp_path, manifest, ("x",))

    assert [r.name for r in results] == ["tests:x"]
    assert results[0].passed
    assert "2" in results[0].detail


def test_a_failing_suite_fails_the_gate_and_is_named(tmp_path: Path) -> None:
    _suite(tmp_path / "test/x", "test-ok.sh", exit_code=0)
    _suite(tmp_path / "test/x", "test-bad.sh", exit_code=1, message="boom")
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())

    results = run_feature_suites(tmp_path, manifest, ("x",))

    assert not results[0].passed
    assert "test-bad.sh" in results[0].detail
    assert "test-ok.sh" not in results[0].detail, "only the failures are worth naming"
    assert "boom" in results[0].detail, "the failure's output is what makes CI diagnosable"


def test_a_dropped_features_suites_are_not_run(tmp_path: Path) -> None:
    """Its scripts are not even in the tree, and crediting them would overstate
    coverage the same way the byte-identity gate refuses to."""
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())

    results = run_feature_suites(tmp_path, manifest, ())

    assert results == []


def test_a_feature_without_a_suite_directory_produces_no_gate(tmp_path: Path) -> None:
    manifest = Manifest(features=(_with_tests("x", None),), withdrawn=())

    assert run_feature_suites(tmp_path, manifest, ("x",)) == []


def test_a_merged_feature_whose_suite_directory_is_missing_fails(tmp_path: Path) -> None:
    """The feature merged, so its scripts should be in the tree. A missing
    directory is a wrong manifest path, not an absence of coverage -- and
    passing silently is exactly how a gate stops meaning anything."""
    manifest = Manifest(features=(_with_tests("x", "test/nope"),), withdrawn=())

    results = run_feature_suites(tmp_path, manifest, ("x",))

    assert not results[0].passed
    assert "test/nope" in results[0].detail


def test_an_empty_suite_directory_fails(tmp_path: Path) -> None:
    """Zero scripts would otherwise report as a pass over nothing."""
    (tmp_path / "test/x").mkdir(parents=True)
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())

    results = run_feature_suites(tmp_path, manifest, ("x",))

    assert not results[0].passed
    assert "no test-*.sh" in results[0].detail


def test_suites_run_with_the_repo_as_the_working_directory(tmp_path: Path) -> None:
    """The suites locate the binary relative to the tree they ship in."""
    d = tmp_path / "test/x"
    d.mkdir(parents=True)
    script = d / "test-cwd.sh"
    script.write_text('#!/bin/sh\ntest "$(pwd)" = "$1" || exit 1\n')
    script.chmod(0o755)
    # The runner passes the repo root as $1, matching test/altlg/'s own convention.
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())

    assert run_feature_suites(tmp_path, manifest, ("x",))[0].passed


def test_run_gates_includes_the_feature_suites(tmp_path: Path) -> None:
    """Wired into run_gates, not left as a function nothing calls."""
    _suite(tmp_path / "test/x", "test-a.sh", exit_code=0)
    manifest = Manifest(features=(_with_tests("x", "test/x"),), withdrawn=())
    stock = stub_aligner(tmp_path / "stock", "@SQ\tSN:chrM\nr1\t0\tchrM\t1\t60\t10M\n")
    cand = stub_aligner(tmp_path / "cand", "@SQ\tSN:chrM\nr1\t0\tchrM\t1\t60\t10M\n")
    fixtures = tmp_path / "fx"
    fixtures.mkdir()
    for f in ("chrM-human.fa.gz", "chrM-read_1.fa.gz", "chrM-read_2.fa.gz"):
        (fixtures / f).write_bytes(b"")

    results = run_gates(
        cand, stock, fixtures, manifest, tmp_path / "wd", merged=("x",), repo=tmp_path
    )

    assert any(r.name == "tests:x" for r in results)
