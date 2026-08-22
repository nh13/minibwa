"""Tests for manifest parsing and validation."""

from pathlib import Path

import pytest

from minibwa_dist.manifest import ManifestError, load_manifest

VALID = """
[[feature]]
name = "inline-appenders"
branch = "feat/inline-appenders"
required = false
output = "identical"
summary = "inline SAM/PAF appenders"
upstream = { pr = 32, status = "rejected" }

[[withdrawn]]
name = "pe-mapq-low-copy-repeats"
branch = "fix/pe-mapq-low-copy-repeats"
report = "reports/2026-07-20-negative-result.md"
summary = "flat AUC; dominated by thresholding"
"""


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "features.toml"
    path.write_text(text)
    return path


def test_loads_features_and_withdrawn(tmp_path: Path) -> None:
    manifest = load_manifest(_write(tmp_path, VALID))
    assert [f.name for f in manifest.features] == ["inline-appenders"]
    assert manifest.features[0].upstream.pr == 32
    assert manifest.features[0].upstream.status == "rejected"
    assert [w.name for w in manifest.withdrawn] == ["pe-mapq-low-copy-repeats"]


def test_rejects_unknown_status(tmp_path: Path) -> None:
    text = VALID.replace('status = "rejected"', 'status = "bogus"')
    with pytest.raises(ManifestError, match="unknown status 'bogus'"):
        load_manifest(_write(tmp_path, text))


def test_rejects_invalid_toml_as_a_manifest_error(tmp_path: Path) -> None:
    """A plain TOML syntax error must surface as `ManifestError`, the one
    exception type every caller already handles as an ordinary operating
    condition (`cli.main` prints one line and exits 1). Left to `tomllib`
    directly, a syntax error raises `tomllib.TOMLDecodeError` instead -- a
    different type that slips past that handling and reaches a caller as an
    unhandled traceback, exactly like an unhandled `OrphanOverPublishedCache`
    did before `cli.main` was taught about it.
    """
    with pytest.raises(ManifestError, match="not valid TOML"):
        load_manifest(_write(tmp_path, "not [valid toml"))


def test_rejects_withdrawn_status_in_features(tmp_path: Path) -> None:
    """A withdrawn feature must never appear in the build. See design section 4.2."""
    text = VALID.replace('status = "rejected"', 'status = "withdrawn"')
    with pytest.raises(ManifestError, match="withdrawn.*must not appear"):
        load_manifest(_write(tmp_path, text))


def test_rejects_graduated_status_in_features(tmp_path: Path) -> None:
    """S2: a feature upstream now carries must not linger in the build.

    Re-merging a graduated feature onto a master that already contains it is a
    near-certain duplicate-application conflict -> optional drop -> a nightly
    issue, forever, until a human edits the file. The manifest must be
    structurally unable to keep building it.
    """
    for status in ("merged", "superseded"):
        text = VALID.replace('status = "rejected"', f'status = "{status}"')
        with pytest.raises(ManifestError, match="upstream now carries this feature"):
            load_manifest(_write(tmp_path, text))


def test_rejects_a_non_boolean_required(tmp_path: Path) -> None:
    """`bool("false")` is True and `bool("")` is False, so a quoting typo in the
    one field only ever written by hand would load with the flag inverted -- and
    a build-gating feature would silently become optional.
    """
    for value in ('"false"', '""', "0"):
        text = VALID.replace("required = false", f"required = {value}")
        with pytest.raises(ManifestError, match="'required' must be a boolean"):
            load_manifest(_write(tmp_path, text))


def test_upstreamable_defaults_true_and_must_be_a_boolean(tmp_path: Path) -> None:
    manifest = load_manifest(_write(tmp_path, VALID))
    assert manifest.features[0].upstreamable is True

    text = VALID.replace("required = false", 'required = false\nupstreamable = "no"')
    with pytest.raises(ManifestError, match="'upstreamable' must be a boolean"):
        load_manifest(_write(tmp_path, text))

    text = VALID.replace("required = false", "required = false\nupstreamable = false")
    assert load_manifest(_write(tmp_path, text)).features[0].upstreamable is False


def test_conditional_requires_condition_and_negative(tmp_path: Path) -> None:
    text = VALID.replace('output = "identical"', 'output = "conditional"')
    with pytest.raises(ManifestError, match="requires 'condition' and 'negative'"):
        load_manifest(_write(tmp_path, text))


def test_rejects_duplicate_names(tmp_path: Path) -> None:
    with pytest.raises(ManifestError, match="duplicate feature name"):
        load_manifest(_write(tmp_path, VALID + VALID.split("[[withdrawn]]")[0]))


def test_withdrawn_requires_report(tmp_path: Path) -> None:
    text = VALID.replace('report = "reports/2026-07-20-negative-result.md"\n', "")
    with pytest.raises(ManifestError, match="missing 'report'"):
        load_manifest(_write(tmp_path, text))


def test_real_manifest_is_valid() -> None:
    """The shipped manifest must always parse, with the ordering assembly depends on.

    Deliberately no feature count: deleting a graduated feature's block is the
    documented, bot-authored routine maintenance action, so a count assertion
    goes red for exactly the change it is supposed to permit. What must hold is
    the structure -- tooling first and required, hub last.
    """
    manifest = load_manifest(Path(__file__).parents[1] / "features.toml")
    assert manifest.features[0].name == "ops-distro", "ops/distro must merge first"
    assert manifest.features[0].required is True
    assert manifest.features[0].upstreamable is False, "the tooling is never offered upstream"
    assert manifest.features[-1].name == "alt-liftgroup", "hub feature must merge last"
    assert manifest.features[-1].required is False
    assert all(f.required is False for f in manifest.features[1:]), (
        "only the tooling entry may abort a build"
    )


def test_a_feature_can_declare_a_test_suite_directory(tmp_path: Path) -> None:
    manifest = load_manifest(
        _write(
            tmp_path,
            """
[[feature]]
name = "a"
branch = "feat/a"
required = false
output = "identical"
tests = "test/altlg"
summary = "a"
upstream = { status = "unsubmitted" }
""",
        )
    )

    assert manifest.features[0].tests == "test/altlg"


def test_the_test_suite_directory_is_optional(tmp_path: Path) -> None:
    """Most features ship no suite of their own; that must stay legal."""
    manifest = load_manifest(
        _write(
            tmp_path,
            """
[[feature]]
name = "a"
branch = "feat/a"
required = false
output = "identical"
summary = "a"
upstream = { status = "unsubmitted" }
""",
        )
    )

    assert manifest.features[0].tests is None
