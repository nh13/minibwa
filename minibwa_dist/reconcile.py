"""Keep `upstream.status` honest, and act on the transitions that matter.

A manifest that can lie is worse than no manifest. The reconciler compares each
feature's recorded status against the live PR state and reports the drift; the
workflow turns that into a PR editing `features.toml`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from minibwa_dist.manifest import GRADUATED_STATUSES, Manifest

# GitHub PR state -> the manifest status it implies.
_PR_STATE_TO_STATUS = {"open": "open", "merged": "merged", "closed": "rejected"}


@dataclass(frozen=True)
class StatusChange:
    """A recorded status that no longer matches reality."""

    feature: str
    old: str
    new: str


def reconcile(manifest: Manifest, pr_states: dict[int, str]) -> list[StatusChange]:
    """Compare recorded statuses against live PR states.

    `pr_states` maps upstream PR number to one of "open", "merged", "closed".
    Features with no PR, or whose PR is absent from the mapping, are left alone.

    There is deliberately no guard preserving a recorded `superseded`: a
    `[[feature]]` can never hold one, because `load_manifest` refuses graduated
    statuses outright and every production caller builds its `Manifest` that way.
    A branch that cannot be reached is not a safety net; it is a claim about the
    data model that stopped being true.
    """
    changes: list[StatusChange] = []
    for feature in manifest.features:
        pr = feature.upstream.pr
        if pr is None or pr not in pr_states:
            continue
        implied = _PR_STATE_TO_STATUS.get(pr_states[pr])
        if implied is None or implied == feature.upstream.status:
            continue
        changes.append(StatusChange(feature=feature.name, old=feature.upstream.status, new=implied))
    return changes


def graduated(manifest: Manifest, changes: list[StatusChange]) -> list[str]:
    """Features that upstream now carries, so they must leave the manifest."""
    return [c.feature for c in changes if c.new in GRADUATED_STATUSES]


def nomination_candidates(manifest: Manifest) -> list[str]:
    """Features never offered upstream and worth offering.

    Withdrawn work is structurally excluded: it does not live in `features`, so
    it can never be nominated. That is the whole reason the status exists.
    Anything the manifest marks `upstreamable = false` is excluded too -- the
    distribution's own tooling is not an upstream contribution. That used to be
    a literal `!= "ops-distro"` here, which quietly re-enabled nominations for
    the tooling the moment it was renamed, and could never cover a second
    never-upstreamable entry.
    """
    return [
        f.name for f in manifest.features if f.upstream.status == "unsubmitted" and f.upstreamable
    ]


_BLOCK_HEADER = re.compile(r"^\[\[(?:feature|withdrawn)\]\]", re.M)
# The lookbehind keeps a key that merely *ends* in the word -- `sub_status`,
# `upstream.status` -- from passing as the bare field.
_STATUS_FIELD = re.compile(r'(?<![\w.\-])(status\s*=\s*")([^"]*)(")')
_NAME_FIELD = r'(?<![\w.\-])name\s*=\s*"{name}"'


def _string_span(text: str, start: int) -> tuple[int, int]:
    """Span of the TOML string literal opening at `start`, delimiters included."""
    quote = text[start]
    delimiter = quote * 3 if text.startswith(quote * 3, start) else quote
    # Only basic (double-quoted) strings honour backslash escapes; in a literal
    # string a backslash is just a backslash.
    escaped = quote == '"'
    index = start + len(delimiter)
    while index < len(text):
        if escaped and text[index] == "\\":
            index += 2
            continue
        if text.startswith(delimiter, index):
            return (start, index + len(delimiter))
        index += 1
    return (start, len(text))  # unterminated: swallow the remainder, never loop


def _inert_spans(text: str) -> list[tuple[int, int]]:
    """Spans of `text` that are string values or comments -- prose, not structure.

    `[[feature]]` and `status = "..."` are ordinary characters inside a quoted
    value, and this manifest carries multi-line `summary` fields written for
    humans -- prose that can quote anything, including the shape of a TOML
    field. A leftmost text search that does not know where strings start and end
    will rewrite `status = "open"` inside a summary and leave the real field
    alone: no exception, no non-zero exit, no diff-shape anomaly, and a manifest
    that quietly keeps lying. Skipping these spans is what makes the leftmost
    match trustworthy without a TOML round-trip library -- `tomllib` reads but
    cannot write, and a writer would be a dependency this project does not take.
    """
    spans: list[tuple[int, int]] = []
    index = 0
    while index < len(text):
        char = text[index]
        if char == "#":
            end = text.find("\n", index)
            spans.append((index, len(text) if end < 0 else end))
        elif char in "\"'":
            spans.append(_string_span(text, index))
        else:
            index += 1
            continue
        index = spans[-1][1]
    return spans


def _find_live(text: str, pattern: re.Pattern[str], start: int, end: int) -> re.Match[str] | None:
    """First match of `pattern` within `text[start:end]` that is real structure.

    Offsets in the returned match are absolute, into `text`.
    """
    inert = _inert_spans(text)
    return next(
        (
            match
            for match in pattern.finditer(text, start, end)
            if not any(lo <= match.start() < hi for lo, hi in inert)
        ),
        None,
    )


def _block_spans(text: str) -> list[tuple[int, int]]:
    """Character spans of each [[feature]] / [[withdrawn]] block."""
    inert = _inert_spans(text)
    starts = [
        match.start()
        for match in _BLOCK_HEADER.finditer(text)
        if not any(lo <= match.start() < hi for lo, hi in inert)
    ]
    bounds = [*starts, len(text)]
    return [(bounds[i], bounds[i + 1]) for i in range(len(starts))]


def _removal_span(text: str, start: int, end: int) -> tuple[int, int]:
    """Narrow a block span to the block itself, so deleting it spares its neighbours.

    `_block_spans` runs each block to the start of the next header, which sweeps
    up the blank line and the section comment that introduce the *next* block --
    `features.toml` has three such headings. Deleting a graduated feature must
    not take an unrelated section comment with it.
    """
    lines = text[start:end].splitlines(keepends=True)
    while lines and (not lines[-1].strip() or lines[-1].lstrip().startswith("#")):
        lines.pop()
    cut = start + sum(len(line) for line in lines)
    # A block sitting between two blank lines would leave both behind. Take one.
    if text[:start].endswith("\n\n") and text[cut:].startswith("\n"):
        cut += 1
    return start, cut


def apply_changes(manifest_text: str, changes: list[StatusChange]) -> str:
    """Return `features.toml` text with each change applied.

    Two shapes, because two shapes are what the manifest will accept. A status
    that merely moved (`open` -> `rejected`) is rewritten in place. A
    *graduation* -- `merged` or `superseded`, meaning upstream now carries the
    feature -- deletes the block outright: design section 7C has the bot open a
    PR deleting the feature, and `load_manifest` rejects those statuses inside
    `[[feature]]` precisely so the build cannot keep re-merging something master
    already has. Writing the status instead would emit a manifest nothing can
    load, and the reconcile workflow would open a PR whose merge wedges sync,
    reconcile and test alike. The record of the graduation belongs in that PR,
    which is where a human reads it -- not in a file the engine must refuse.

    The invariant, enforced by test: whatever this returns must load.

    Operates on raw text, not a parsed document, so comments and formatting
    survive -- the manifest is read by humans as much as by the engine.

    Edits are confined to the target feature's own block. A regex spanning
    `name = "<x>" ... status = "<old>"` with re.S looks equivalent but is not:
    when `old` does not match what the file actually holds, the scan runs past
    the block and rewrites the NEXT feature that happens to carry `old` -- one
    substitution, so a count-based guard never fires, and an unrelated feature
    is silently corrupted. Locate the block, verify, then act.
    """
    text = manifest_text
    for change in changes:
        name_field = re.compile(_NAME_FIELD.format(name=re.escape(change.feature)))
        target = next(
            (span for span in _block_spans(text) if _find_live(text, name_field, *span)),
            None,
        )
        if target is None:
            raise ValueError(
                f"no block declares feature '{change.feature}'; edit features.toml by hand"
            )

        start, end = target
        found = _find_live(text, _STATUS_FIELD, start, end)
        if found is None or found.group(2) != change.old:
            actual = found.group(2) if found else "<no status field>"
            raise ValueError(
                f"'{change.feature}' has status '{actual}', not '{change.old}' "
                f"(wanted {change.old} -> {change.new}); edit features.toml by hand"
            )

        if change.new in GRADUATED_STATUSES:
            cut_start, cut_end = _removal_span(text, start, end)
            text = text[:cut_start] + text[cut_end:]
        else:
            text = text[: found.start(2)] + change.new + text[found.end(2) :]
    return text
