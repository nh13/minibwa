"""Tests for label reconciliation logic."""

from pathlib import Path

import pytest

from minibwa_dist.manifest import Feature, Manifest, Upstream, Withdrawn, load_manifest
from minibwa_dist.reconcile import (
    StatusChange,
    apply_changes,
    graduated,
    nomination_candidates,
    reconcile,
)

_MANIFEST = Path(__file__).resolve().parents[1] / "features.toml"

# A stand-in for the shipped manifest's shape: section comments introducing
# blocks, a feature with no status field at all, and later blocks that do carry
# one. Inline rather than read off `minibwa_dist/features.toml`, which the reconcile
# bot in this same repo rewrites on a schedule -- a test that asserts on a live,
# machine-edited file goes red for reasons that have nothing to do with the code
# under test, and the bot's own PR is what turns it red.
_SYNTHETIC = """\
# --- a heading that introduces the block below ---
[[feature]]
name = "no-status"
branch = "feat/no-status"
required = false
output = "identical"
summary = "never submitted, so it carries no status field of its own"
upstream = { status = "unsubmitted" }

# --- a heading that belongs to the block below it ---
[[feature]]
name = "tracked"
branch = "feat/tracked"
required = false
output = "identical"
summary = "an open PR upstream"
upstream = { pr = 20, status = "open" }

# --- a heading that must survive a deletion of the block above it ---
[[feature]]
name = "also-tracked"
branch = "feat/also-tracked"
required = false
output = "identical"
summary = "another open PR upstream"
upstream = { pr = 21, status = "open" }
"""


def _load_text(tmp_path: Path, text: str) -> Manifest:
    """Parse `text` as a manifest, the way the engine would after a bot rewrote it."""
    path = tmp_path / "features.toml"
    path.write_text(text)
    return load_manifest(path)


def _feature(
    name: str, status: str, pr: int | None = None, *, upstreamable: bool = True
) -> Feature:
    return Feature(
        name=name,
        branch=f"feat/{name}",
        required=False,
        output="identical",
        summary=name,
        upstream=Upstream(status=status, pr=pr),
        upstreamable=upstreamable,
    )


def test_detects_a_merged_pr() -> None:
    manifest = Manifest(features=(_feature("a", "open", 33),), withdrawn=())
    changes = reconcile(manifest, {33: "merged"})
    assert changes == [StatusChange(feature="a", old="open", new="merged")]
    assert graduated(manifest, changes) == ["a"]


def test_detects_a_closed_pr() -> None:
    manifest = Manifest(features=(_feature("a", "open", 33),), withdrawn=())
    changes = reconcile(manifest, {33: "closed"})
    assert changes[0].new == "rejected"
    assert graduated(manifest, changes) == [], "a rejection is not a graduation"


def test_no_change_when_state_already_matches() -> None:
    manifest = Manifest(features=(_feature("a", "rejected", 32),), withdrawn=())
    assert reconcile(manifest, {32: "closed"}) == []


def test_features_without_a_pr_are_untouched() -> None:
    manifest = Manifest(features=(_feature("a", "unsubmitted"),), withdrawn=())
    assert reconcile(manifest, {}) == []


def test_nominations_exclude_withdrawn_work_and_tooling() -> None:
    """The bug `withdrawn` exists to prevent: never re-offer disproven work.

    The tooling is excluded by the manifest's own `upstreamable = false`, not by
    its name: keyed on the name, renaming the entry silently started nominating
    the distribution's own machinery to upstream.
    """
    manifest = Manifest(
        features=(
            _feature("the-tooling-whatever-it-is-called", "unsubmitted", upstreamable=False),
            _feature("fresh", "unsubmitted"),
            _feature("done", "rejected"),
        ),
        withdrawn=(Withdrawn(name="dead", branch="fix/dead", report="r.md", summary="disproven"),),
    )
    assert nomination_candidates(manifest) == ["fresh"]


def test_apply_changes_rewrites_status_and_keeps_comments() -> None:
    text = """
# a comment that must survive
[[feature]]
name = "a"
branch = "feat/a"
upstream = { pr = 33, status = "open" }
"""
    out = apply_changes(text, [StatusChange(feature="a", old="open", new="rejected")])
    assert 'status = "rejected"' in out
    assert "# a comment that must survive" in out


def test_apply_changes_refuses_a_miss() -> None:
    """A silent no-op would let the manifest keep lying."""
    with pytest.raises(ValueError, match="has status"):
        apply_changes(
            '[[feature]]\nname = "a"\n', [StatusChange(feature="a", old="open", new="merged")]
        )


def test_apply_changes_does_not_bleed_into_the_next_block() -> None:
    """Regression: a non-greedy `.*?` scan with no block boundary would run past
    a feature with no matching status and silently rewrite the *next* feature
    downstream that happens to carry `old`. `no-status` has no `status = "open"`
    anywhere in its block; `tracked`, two blocks later, does. Asking to change
    `no-status` from `old="open"` must raise, not quietly corrupt `tracked`.
    """
    with pytest.raises(ValueError, match="has status"):
        apply_changes(
            _SYNTHETIC,
            [StatusChange(feature="no-status", old="open", new="merged")],
        )

    tracked_block = _SYNTHETIC.split('name = "tracked"', 1)[1]
    assert 'status = "open"' in tracked_block.split("[[feature]]", 1)[0]


def test_apply_changes_rejects_an_unknown_feature() -> None:
    with pytest.raises(ValueError, match="no block declares"):
        apply_changes(
            '[[feature]]\nname = "a"\nstatus = "open"\n',
            [StatusChange(feature="ghost", old="open", new="merged")],
        )


def test_apply_changes_picks_the_right_block_among_shared_statuses() -> None:
    text = '[[feature]]\nname = "a"\nstatus = "open"\n\n[[feature]]\nname = "b"\nstatus = "open"\n'
    out = apply_changes(text, [StatusChange(feature="b", old="open", new="rejected")])

    a_block = out.split('name = "a"', 1)[1].split("[[feature]]", 1)[0]
    b_block = out.split('name = "b"', 1)[1]
    assert 'status = "open"' in a_block
    assert 'status = "rejected"' in b_block


def test_apply_changes_ignores_a_status_shaped_string_in_the_prose() -> None:
    """Regression: `summary` is free-form prose and TOML literal strings need no
    escaping, so a summary can legitimately contain the exact text of a field.
    A leftmost text search rewrote *that* and left `upstream.status` alone --
    no exception, no non-zero exit, and a manifest still lying about the
    feature while an unrelated sentence had been silently edited.
    """
    text = (
        "[[feature]]\n"
        'name = "quoter"\n'
        "summary = 'their API used to report status = \"open\" for this PR'\n"
        'upstream = { pr = 41, status = "open" }\n'
    )
    out = apply_changes(text, [StatusChange(feature="quoter", old="open", new="rejected")])

    assert 'their API used to report status = "open" for this PR' in out, "prose was rewritten"
    assert 'upstream = { pr = 41, status = "rejected" }' in out


def test_apply_changes_ignores_a_status_shaped_string_in_a_multiline_summary() -> None:
    """The real manifest's `[[withdrawn]]` entries carry multi-line `summary`
    values, so the string scanner has to survive `\"\"\"` blocks -- including a
    line inside one that starts like a block header.
    """
    text = (
        "[[feature]]\n"
        'name = "quoter"\n'
        'summary = """we tried\n'
        '[[feature]] status = "merged"\n'
        'was the note upstream left"""\n'
        'upstream = { pr = 41, status = "open" }\n'
    )
    out = apply_changes(text, [StatusChange(feature="quoter", old="open", new="rejected")])

    assert '[[feature]] status = "merged"' in out, "prose was rewritten"
    assert 'upstream = { pr = 41, status = "rejected" }' in out


def test_graduation_deletes_the_block_and_the_result_still_loads(tmp_path: Path) -> None:
    """The invariant: what `apply_changes` writes, `load_manifest` must read.

    Regression: writing `status = "merged"` produced a manifest `load_manifest`
    refuses outright, so the reconcile bot would open a PR whose merge wedges
    sync, reconcile and test alike. Design section 7C says delete the block.
    """
    before = _load_text(tmp_path, _SYNTHETIC)

    out = apply_changes(_SYNTHETIC, [StatusChange(feature="tracked", old="open", new="merged")])
    after = _load_text(tmp_path, out)

    assert {f.name for f in after.features} == {f.name for f in before.features} - {"tracked"}
    assert after.withdrawn == before.withdrawn
    # A section heading introduces the block *below* it, so deleting a block must
    # leave the next block's heading -- and every other -- exactly where it was.
    assert "# --- a heading that must survive a deletion of the block above it ---" in out
    assert "# --- a heading that introduces the block below ---" in out
    assert 'name = "tracked"' not in out


def test_every_reconcilable_change_produces_a_loadable_manifest(tmp_path: Path) -> None:
    """Exhaustive over what `reconcile()` can emit.

    `reconcile` derives changes from PR state alone, so every feature carrying a
    PR can be pushed to any of the three implied statuses. All of them must
    round-trip through `load_manifest` -- a graduation by deleting the block, a
    move by rewriting the field.
    """
    manifest = _load_text(tmp_path, _SYNTHETIC)
    tracked = [f for f in manifest.features if f.upstream.pr is not None]
    assert len(tracked) == 2, "the fixture must keep covering more than one tracked block"

    checked = 0
    for feature in tracked:
        for new in ("open", "merged", "rejected"):
            if new == feature.upstream.status:
                continue
            change = StatusChange(feature=feature.name, old=feature.upstream.status, new=new)
            _load_text(tmp_path, apply_changes(_SYNTHETIC, [change]))
            checked += 1
    assert checked == 4


def test_the_shipped_manifest_round_trips_through_reconcile(tmp_path: Path) -> None:
    """A smoke test against the real file, asserting only what stays true as the
    bot edits it: whatever `reconcile` could emit for it must still load.
    """
    original = _MANIFEST.read_text()
    for feature in load_manifest(_MANIFEST).features:
        if feature.upstream.pr is None:
            continue
        for new in ("open", "merged", "rejected"):
            if new == feature.upstream.status:
                continue
            change = StatusChange(feature=feature.name, old=feature.upstream.status, new=new)
            _load_text(tmp_path, apply_changes(original, [change]))
