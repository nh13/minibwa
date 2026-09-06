"""Generate the distribution's documentation from the manifest.

The manifest is the source of truth; every human-readable list of what this
fork contains is derived from it, so the two cannot drift.
"""

from __future__ import annotations

from pathlib import Path

from minibwa_dist.manifest import Manifest

UPSTREAM = "lh3/minibwa"
BLOCK_BEGIN = "<!-- DISTRO:BEGIN -->"
BLOCK_END = "<!-- DISTRO:END -->"

# The benchmark CTA is generated OUT OF TREE by the bench harness
# (minibwa-ab-bench: `python -m cta.bench render`) and committed here as a
# self-contained artifact. This file consumes it; it never measures anything.
# Keep `render_cta` in sync with that harness's cta/render_cta.py -- the artifact
# format below is the contract between them. See CLAUDE.md "The CTA design".
CTA_ARTEFACT = "minibwa_dist/cta.md"
CTA_FEATURES_MARKER = "<!-- CTA:FEATURES "
CTA_DATE_MARKER = "<!-- CTA:MEASURED "
CTA_MARKER_END = " -->"

BANNER = """\
> **This is a downstream build of [lh3/minibwa](https://github.com/lh3/minibwa).**
> It carries changes upstream declined, plus a few not yet offered. Upstream remains the source
> of truth for everything else.
>
> **Install from a release tag, not from a branch.** Tags such as `v0.6-nh13.1` are immutable;
> the `dist` branch is rebuilt and force-pushed on every upstream change and will rewind under
> you. `minibwa version` reports the downstream version, and it appears in the `@PG VN:` tag of
> every SAM file this build writes, so output is always traceable to the build that produced it.
>
> Bug reports for anything in the table below belong here, not upstream.
> Changes investigated and deliberately not shipped are in [`GRAVEYARD.md`](GRAVEYARD.md).
"""


def _pr_ref(pr: int | None) -> str:
    return f"[{UPSTREAM}#{pr}](https://github.com/{UPSTREAM}/pull/{pr})" if pr else "—"


def render_cta(cta_path: Path | None, manifest: Manifest) -> str:
    """The benchmark headline for the README front page, or "" / a pending line.

    The artifact records the feature set the numbers were measured against. If
    that set matches the manifest being built, the numbers are shown; if it
    drifted (a feature added or removed since), they are hidden behind a
    one-line "pending re-measurement" note -- the hard-hide rule, so the front
    page never advertises numbers for a build it no longer describes. A missing
    artifact (never measured) renders nothing.

    A names-level check: a same-name feature re-cut is not caught here; the
    measured date in the artifact's caption carries that staleness instead.

    The result is a pure function of (artifact bytes, manifest feature names) --
    no SHA, date, or clock -- so a nightly re-render is byte-stable and does not
    trip the sync's no-op detector.
    """
    if cta_path is None or not cta_path.exists():
        return ""
    measured: set[str] = set()
    date = "an earlier run"
    body: list[str] = []
    for line in cta_path.read_text().splitlines():
        if line.startswith(CTA_FEATURES_MARKER) and line.endswith(CTA_MARKER_END):
            names = line[len(CTA_FEATURES_MARKER) : -len(CTA_MARKER_END)]
            measured = {n.strip() for n in names.split(",") if n.strip()}
        elif line.startswith(CTA_DATE_MARKER) and line.endswith(CTA_MARKER_END):
            date = line[len(CTA_DATE_MARKER) : -len(CTA_MARKER_END)].strip()
        else:
            body.append(line)
    if measured and measured == {f.name for f in manifest.features}:
        return "\n".join(body).strip()
    return (
        "> ⚡ **Benchmark numbers are hidden pending re-measurement** — the "
        f"distribution's feature set has changed since they were last measured ({date})."
    )


def render_feature_table(manifest: Manifest) -> str:
    """A markdown table of everything the distribution adds to upstream."""
    lines = [
        "| feature | upstream | status | output | summary |",
        "|---|---|---|---|---|",
    ]
    for feature in manifest.features:
        lines.append(
            f"| `{feature.name}` | {_pr_ref(feature.upstream.pr)} "
            f"| {feature.upstream.status} | {feature.output} | {feature.summary} |"
        )
    return "\n".join(lines)


def render_graveyard(manifest: Manifest) -> str:
    """Investigated and deliberately not shipped, each with its disproof."""
    parts = [
        "# Graveyard",
        "",
        "Changes investigated and deliberately **not** shipped. They are recorded so the",
        "same idea is not re-litigated, and their branches are kept so each disproof stays",
        "reproducible: the branch is the artifact you can check out and re-measure.",
        "Nothing here is a candidate for upstreaming.",
        "",
        "The report path on each entry names the maintainer's own analysis notes. Those",
        "are not part of this repository — the branch and this summary are what travel",
        "with the distribution.",
        "",
    ]
    for item in manifest.withdrawn:
        parts += [
            f"## `{item.name}`",
            "",
            f"- **Branch:** `{item.branch}`",
            f"- **Report:** `{item.report}`",
            "",
            item.summary,
            "",
        ]
    return "\n".join(parts)


def render_release_notes(
    manifest: Manifest, upstream_notes: str, upstream_tag: str, dropped: tuple[str, ...]
) -> str:
    """Release notes: upstream's own, plus what this distribution adds and omits."""
    parts = [
        f"Downstream build of [{UPSTREAM} {upstream_tag}]"
        f"(https://github.com/{UPSTREAM}/releases/tag/{upstream_tag}).",
        "",
        "## Upstream release notes",
        "",
        upstream_notes,
        "",
        "## Added by this distribution",
        "",
        render_feature_table(manifest),
        "",
    ]
    if dropped:
        parts += [
            "## Features NOT included in this build",
            "",
            "These are in the manifest but would not merge against this upstream revision.",
            "They will return once the conflict is resolved:",
            "",
            *[f"- `{name}`" for name in dropped],
            "",
        ]
    return "\n".join(parts)


def update_readme(readme: Path, manifest: Manifest, cta_path: Path | None = None) -> None:
    """Insert or refresh the distribution banner and feature table in `readme`.

    The block is delimited so repeated renders replace rather than accumulate, and
    upstream's own prose is preserved around it. `dist` is the default branch, so
    this file is the repo's front page -- it must say what this fork is.

    `cta_path`, when given and present, prepends the benchmark headline above the
    banner (see `render_cta`); when absent or drifted it contributes nothing or a
    pending line, so the block is always well-formed.
    """
    cta = render_cta(cta_path, manifest)
    cta_part = f"{cta}\n\n" if cta else ""
    block = f"{BLOCK_BEGIN}\n{cta_part}{BANNER}\n{render_feature_table(manifest)}\n{BLOCK_END}"
    text = readme.read_text() if readme.exists() else ""

    begins = text.count(BLOCK_BEGIN)
    ends = text.count(BLOCK_END)
    if begins == 1 and ends == 1 and text.index(BLOCK_BEGIN) < text.index(BLOCK_END):
        head, _, rest = text.partition(BLOCK_BEGIN)
        _, _, tail = rest.partition(BLOCK_END)
        readme.write_text(f"{head}{block}{tail}")
        return

    if begins or ends:
        # Anything but exactly one well-ordered pair is corruption -- a truncated
        # file, or a merge that landed half the block or two copies of it. Both
        # of the other shapes used to slip past a `BEGIN in text and END in text`
        # guard into the rewrite: with END before BEGIN, nothing after BEGIN
        # contains END, so `partition` returned an empty tail and every byte
        # after BEGIN was deleted with no diagnostic; with two complete pairs,
        # the second survived every render, contradicting this function's whole
        # reason to delimit. Refuse instead: this file is the front page of the
        # default branch, rewritten unattended on every sync.
        found = (
            "they appear in the wrong order"
            if begins == 1 and ends == 1
            else f"found {begins} begin and {ends} end marker(s)"
        )
        raise ValueError(
            f"{readme}: expected exactly one {BLOCK_BEGIN} ... {BLOCK_END} pair, but "
            f"{found} -- the distribution block is corrupt; repair it by hand "
            "before re-rendering"
        )

    # First insertion: after the title line if there is one, else at the top.
    lines = text.splitlines(keepends=True)
    at = 1 if lines and lines[0].startswith("#") else 0
    readme.write_text("".join(lines[:at]) + f"\n{block}\n\n" + "".join(lines[at:]))
