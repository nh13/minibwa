"""Tests for generated documentation."""

from pathlib import Path

import pytest

from minibwa_dist.manifest import Feature, Manifest, Upstream, Withdrawn
from minibwa_dist.render import (
    BLOCK_BEGIN,
    BLOCK_END,
    render_feature_table,
    render_graveyard,
    render_release_notes,
)

MANIFEST = Manifest(
    features=(
        Feature(
            name="inline-appenders",
            branch="feat/inline-appenders",
            required=False,
            output="identical",
            summary="inline SAM/PAF appenders",
            upstream=Upstream(status="rejected", pr=32),
        ),
    ),
    withdrawn=(
        Withdrawn(
            name="pe-mapq",
            branch="fix/pe-mapq-low-copy-repeats",
            report="reports/negative.md",
            summary="dominated by thresholding",
        ),
    ),
)


def test_feature_table_links_the_upstream_pr() -> None:
    table = render_feature_table(MANIFEST)
    assert "inline-appenders" in table
    assert "lh3/minibwa#32" in table
    assert "rejected" in table


def test_graveyard_names_the_report() -> None:
    text = render_graveyard(MANIFEST)
    assert "reports/negative.md" in text
    assert "pe-mapq" in text


def test_release_notes_list_dropped_features_loudly() -> None:
    """A silent drop is the difference between a distribution and a surprise."""
    notes = render_release_notes(MANIFEST, "upstream notes here", "v0.6", ("alt-liftgroup",))
    assert "upstream notes here" in notes
    assert "alt-liftgroup" in notes
    assert "not included" in notes.lower()


def test_release_notes_omit_the_dropped_section_when_none() -> None:
    notes = render_release_notes(MANIFEST, "notes", "v0.6", ())
    assert "not included" not in notes.lower()


def test_update_readme_replaces_the_block_idempotently(tmp_path: Path) -> None:
    """The front page must regenerate, not drift, when features.toml changes."""
    from minibwa_dist.render import update_readme

    readme = tmp_path / "README.md"
    readme.write_text("# minibwa\n\nUpstream prose.\n")
    update_readme(readme, MANIFEST)
    once = readme.read_text()
    assert "downstream build" in once.lower()
    assert "inline-appenders" in once
    assert "Upstream prose." in once, "upstream content must be preserved"

    update_readme(readme, MANIFEST)
    assert readme.read_text() == once, "re-rendering must not duplicate the block"


def test_update_readme_refuses_a_begin_without_an_end(tmp_path: Path) -> None:
    """A truncated file must not be silently 'fixed' by inserting a second block."""
    from minibwa_dist.render import update_readme

    readme = tmp_path / "README.md"
    original = f"# minibwa\n\n{BLOCK_BEGIN}\nstray content, no closing marker\n"
    readme.write_text(original)

    with pytest.raises(ValueError, match="exactly one"):
        update_readme(readme, MANIFEST)

    assert readme.read_text() == original, "the corrupt file must be left untouched"


def test_update_readme_refuses_an_end_without_a_begin(tmp_path: Path) -> None:
    """Same defect, mirrored: an orphaned END must also be refused, not repaired."""
    from minibwa_dist.render import update_readme

    readme = tmp_path / "README.md"
    original = f"# minibwa\n\nstray content, no opening marker\n{BLOCK_END}\n"
    readme.write_text(original)

    with pytest.raises(ValueError, match="exactly one"):
        update_readme(readme, MANIFEST)

    assert readme.read_text() == original, "the corrupt file must be left untouched"


def test_update_readme_refuses_markers_in_the_wrong_order(tmp_path: Path) -> None:
    """Silent data loss: END before BEGIN satisfied a `both present` guard, but
    nothing after BEGIN then contained END, so `partition` returned an empty tail
    and every byte after BEGIN was deleted -- no exception, no diagnostic, on the
    front page of the default branch, rewritten unattended on every sync.
    """
    from minibwa_dist.render import update_readme

    readme = tmp_path / "README.md"
    original = f"# minibwa\n\nPROSE A\n{BLOCK_END}\nPROSE B\n{BLOCK_BEGIN}\nPROSE C\n"
    readme.write_text(original)

    with pytest.raises(ValueError, match="wrong order"):
        update_readme(readme, MANIFEST)

    assert readme.read_text() == original, "nothing may be dropped from a file we refuse"


def test_update_readme_refuses_a_duplicated_block(tmp_path: Path) -> None:
    """Two complete pairs: only the first was ever replaced, so the stale copy
    survived every render -- the exact accumulation the markers exist to prevent.
    """
    from minibwa_dist.render import update_readme

    readme = tmp_path / "README.md"
    block = f"{BLOCK_BEGIN}\nold\n{BLOCK_END}\n"
    original = f"# minibwa\n\n{block}\nprose\n\n{block}"
    readme.write_text(original)

    with pytest.raises(ValueError, match="2 begin and 2 end"):
        update_readme(readme, MANIFEST)

    assert readme.read_text() == original


def test_graveyard_says_the_report_is_not_in_the_repository() -> None:
    """The paths name the maintainer's own notes, which never travel with `dist`;
    GRAVEYARD.md ships on the public default branch and used to read as if a
    reader could open them.
    """
    text = render_graveyard(MANIFEST)
    assert "not part of this repository" in text
