"""Parse and validate `features.toml`, the distribution manifest.

The manifest is the single source of truth for what the distribution contains.
Validation is strict on purpose: a manifest that parses is a manifest the
assembly engine can execute without further checks.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

VALID_STATUSES = frozenset({"unsubmitted", "open", "rejected", "merged", "superseded", "withdrawn"})
VALID_OUTPUTS = frozenset({"identical", "conditional", "changes-output"})

# Statuses that mean "upstream now has this", so the feature must leave the build.
GRADUATED_STATUSES = frozenset({"merged", "superseded"})


class ManifestError(ValueError):
    """Raised when the manifest is malformed or internally inconsistent."""


@dataclass(frozen=True)
class Upstream:
    """Upstream tracking state for a feature."""

    status: str
    pr: int | None = None


@dataclass(frozen=True)
class Feature:
    """One feature that composes the distribution."""

    name: str
    branch: str
    required: bool
    output: str
    summary: str
    upstream: Upstream
    condition: str | None = None
    negative: str | None = None
    # Directory of `test-*.sh` suites this feature ships, run by the output gates
    # when -- and only when -- the feature is in the build. This is the positive-
    # path coverage the byte-identity gate deliberately does not provide: those
    # paths change output on purpose, so identity says nothing about them.
    tests: str | None = None
    # False for work that is ours to carry and never upstream's to take -- the
    # distribution's own tooling. Defaults true: a feature is a candidate for
    # upstreaming unless the manifest says otherwise.
    upstreamable: bool = True


@dataclass(frozen=True)
class Withdrawn:
    """A feature investigated and deliberately not shipped. Never built."""

    name: str
    branch: str
    report: str
    summary: str


@dataclass(frozen=True)
class Manifest:
    """The parsed manifest."""

    features: tuple[Feature, ...]
    withdrawn: tuple[Withdrawn, ...]


def _require(table: dict[str, object], key: str, where: str) -> object:
    if key not in table:
        raise ManifestError(f"{where}: missing '{key}'")
    return table[key]


def _require_bool(table: dict[str, object], key: str, where: str) -> bool:
    """A TOML boolean, refusing anything `bool()` would coerce.

    `bool("false")` is True and `bool("")` is False, so a quoting typo --
    `required = "false"`, the one field only ever written by hand -- would load
    cleanly with the flag inverted, and a build-gating feature would silently
    become optional. Every other field is checked against a set of allowed
    values; this one was not checked at all.
    """
    value = _require(table, key, where)
    if not isinstance(value, bool):
        raise ManifestError(f"{where}: '{key}' must be a boolean, not {value!r}")
    return value


def _optional_bool(table: dict[str, object], key: str, where: str, default: bool) -> bool:
    return _require_bool(table, key, where) if key in table else default


def _parse_feature(table: dict[str, object], index: int) -> Feature:
    where = f"feature[{index}]"
    name = str(_require(table, "name", where))
    where = f"feature '{name}'"

    upstream_raw = _require(table, "upstream", where)
    if not isinstance(upstream_raw, dict):
        raise ManifestError(f"{where}: 'upstream' must be a table")
    status = str(_require(upstream_raw, "status", where))
    if status not in VALID_STATUSES:
        raise ManifestError(f"{where}: unknown status '{status}'")
    if status == "withdrawn":
        raise ManifestError(
            f"{where}: withdrawn features must not appear in [[feature]]; "
            "use [[withdrawn]] so the build never sees them"
        )
    if status in GRADUATED_STATUSES:
        # Upstream taking a feature means master already contains it, so
        # re-merging it every sync is a near-certain duplicate-application
        # conflict -> optional drop -> a nightly issue, forever, until a human
        # notices and edits the file by hand. `reconcile` already reports a
        # graduation; this makes the manifest structurally unable to keep
        # building a feature upstream now carries.
        raise ManifestError(
            f"{where}: status '{status}' means upstream now carries this feature; "
            "delete its [[feature]] block instead of leaving it in the manifest"
        )

    output = str(_require(table, "output", where))
    if output not in VALID_OUTPUTS:
        raise ManifestError(f"{where}: unknown output class '{output}'")

    condition = table.get("condition")
    negative = table.get("negative")
    tests = table.get("tests")
    if output == "conditional" and not (condition and negative):
        # Documentation, and required as such: the gate runs one default-flags
        # comparison that stands in for every conditional feature's negative
        # case at once, so nothing reads these strings. They record which
        # negative case that single run is standing in for, which is the only
        # way a reader can check the gate actually covers the feature.
        raise ManifestError(
            f"{where}: output='conditional' requires 'condition' and 'negative' "
            "so the negative case the byte-identity gate stands in for is recorded"
        )

    pr_raw = upstream_raw.get("pr")
    return Feature(
        name=name,
        branch=str(_require(table, "branch", where)),
        required=_require_bool(table, "required", where),
        output=output,
        summary=str(_require(table, "summary", where)),
        upstream=Upstream(status=status, pr=int(pr_raw) if pr_raw is not None else None),
        condition=str(condition) if condition else None,
        negative=str(negative) if negative else None,
        tests=str(tests) if tests else None,
        upstreamable=_optional_bool(table, "upstreamable", where, default=True),
    )


def _parse_withdrawn(table: dict[str, object], index: int) -> Withdrawn:
    where = f"withdrawn[{index}]"
    name = str(_require(table, "name", where))
    where = f"withdrawn '{name}'"
    return Withdrawn(
        name=name,
        branch=str(_require(table, "branch", where)),
        report=str(_require(table, "report", where)),
        summary=str(_require(table, "summary", where)),
    )


def parse_manifest(text: str) -> Manifest:
    """Parse and validate manifest TOML text, without requiring it be on disk.

    `load_manifest` is `parse_manifest(path.read_text())` plus the file read;
    this half is the one a caller needs to validate text it is *about* to
    write, before it writes it -- e.g. `reconcile --write`, which must refuse
    to put a manifest on disk that this loader cannot read back.

    Raises `ManifestError` on any structural or semantic problem -- including a
    plain TOML syntax error, which `tomllib` itself would otherwise raise as
    `tomllib.TOMLDecodeError`. Every caller of this module already handles
    `ManifestError` as an ordinary operating condition (`cli.main` prints one
    line and exits 1 rather than a traceback); a `TOMLDecodeError` slipping
    past that as a different exception type defeats it just as surely as not
    handling anything at all.
    """
    try:
        raw = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ManifestError(f"not valid TOML: {exc}") from exc

    features = tuple(_parse_feature(t, i) for i, t in enumerate(raw.get("feature", [])))
    withdrawn = tuple(_parse_withdrawn(t, i) for i, t in enumerate(raw.get("withdrawn", [])))

    seen: set[str] = set()
    for item in (*features, *withdrawn):
        if item.name in seen:
            raise ManifestError(f"duplicate feature name '{item.name}'")
        seen.add(item.name)

    return Manifest(features=features, withdrawn=withdrawn)


def load_manifest(path: Path) -> Manifest:
    """Load and validate the manifest at `path`.

    Raises `ManifestError` on any structural or semantic problem.
    """
    return parse_manifest(path.read_text())
