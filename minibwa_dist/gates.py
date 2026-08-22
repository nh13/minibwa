"""Output gates: prove the distribution still aligns exactly like upstream.

Every `identical` feature must leave the SAM byte-for-byte unchanged, and every
`conditional` feature must do so in its negative case (the condition unmet).
This is what makes a rerere-replayed merge safe to trust: rerere reapplies text,
not meaning, and only a real alignment comparison can catch a resolution that is
textually plausible and semantically wrong.
"""

from __future__ import annotations

import hashlib
import subprocess
from dataclasses import dataclass
from pathlib import Path

from minibwa_dist.manifest import Manifest

# @PG carries MB_VERSION and the command line, and @CO is free-form; both differ
# between a stock and a downstream build by construction. Every other header
# line -- @HD, @SQ, @RG -- is real output and a change in it is a regression.
_EXEMPT_HEADERS = (b"@PG", b"@CO")

# One entry per gated invocation: a label, the flags added after `map`, and
# whether the second read file is passed. `default-flags` comes first and is the
# invocation `default-flags-byte-identity` is named for -- it is also the only one
# whose detail carries the feature coverage list, so the other three do not repeat
# a fourteen-name list in the CI log.
#
# The extra four exist because one default paired run leaves most of `format.c`
# unexercised, and `format.c` is the file downstream churns hardest: resolving
# `inline-appenders` required rewriting upstream's new mate-unmapped `r_pri` logic
# in appender style rather than picking a side, which is precisely the shape of
# resolution rerere can replay as textually plausible and semantically wrong. `-u`
# reaches unmapped-record suppression, `-b MD` reaches tag emission, `--eqx`
# reaches CIGAR construction, and the single-end run reaches the unpaired paths a
# paired-only comparison never enters. All are upstream flags with no downstream
# feature enabled, so the `identical` contract binds on every one of them.
#
# Every entry is verified to CHANGE the output on upstream's chrM fixture -- a
# mode whose SAM matches the default one is a second copy of the default gate
# wearing a different name. Two candidates were rejected on exactly that test:
# `-a`, which in `map` is only `flag &= ~MB_F_PAF` and so is a no-op because SAM
# is already the default (it means all-hits in the OTHER subcommand, which is
# where the habit comes from), and `--outn`/`--outs`/`-N`, because this fixture
# produces no secondary alignments for them to emit. Secondary-record emission is
# therefore NOT covered here, and cannot be without a new fixture.
_MODES: tuple[tuple[str, tuple[str, ...], bool], ...] = (
    ("default-flags", (), True),
    ("no-unmapped", ("-u",), True),
    ("base-tag", ("-b", "MD"), True),
    ("eqx-cigar", ("--eqx",), True),
    ("single-end", (), False),
)


# Bounded by the job, not by generosity: both gate jobs in distro-sync.yml cap at
# `timeout-minutes: 60`, and this budget applies to each of a feature's suites in
# turn. A value the whole set can outlast is a guard that never fires -- the
# runner kills the job first and prints nothing at all. 120s x 20 suites is 40
# minutes, leaving the builds and the byte-identity comparison the rest; the
# in-tree ALT suites finish in seconds, so this is still enormous slack per
# suite. `test_the_suite_timeout_fits_inside_the_job_timeout` holds the bound.
_SUITE_TIMEOUT_S = 120

# Enough failures to see the pattern, not so many that one gate line is a log
# dump: 15 failing suites at 800 characters each is a 12KB single line.
_MAX_REPORTED_FAILURES = 3


def _tail(text: str, limit: int = 800) -> str:
    """The last `limit` characters, which is where a shell suite says why it failed."""
    text = text.strip()
    return text if len(text) <= limit else "..." + text[-limit:]


def _gate_name(label: str) -> str:
    """Gate name for a mode label.

    `default-flags-byte-identity` is load-bearing text: the runbook and the
    release notes name it, so it keeps its spelling rather than becoming
    `byte-identity:default-flags` for symmetry.
    """
    return "default-flags-byte-identity" if label == "default-flags" else f"byte-identity:{label}"


@dataclass(frozen=True)
class GateResult:
    """Outcome of one gate."""

    name: str
    passed: bool
    detail: str


def sam_digest(path: Path) -> str:
    """sha256 over a SAM file, excluding only the @PG and @CO header lines."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            if not line.startswith(_EXEMPT_HEADERS):
                digest.update(line)
    return digest.hexdigest()


def align_fixture(binary: Path, fixture_dir: Path, workdir: Path) -> dict[str, Path]:
    """Index the chrM fixture with `binary` and align it once per `_MODES` entry.

    Returns the SAM path for each mode label. One index serves every mode: the
    index is a function of the reference and the binary, not of the mapping flags.

    `minibwa map` takes an index, not a FASTA, so the index is built here -- with
    the SAME binary under test, never shared between stock and candidate, so an
    index-format regression from `index-threads-v2` cannot hide behind a stock index.

    Uses upstream's own `test/chrM-*.fa.gz`; no fixture data is added to the repo.
    """
    # `Path("./minibwa")` normalises to `"minibwa"` (pathlib drops the leading
    # `./`), and subprocess then treats a slash-free string as a PATH lookup
    # rather than a path relative to cwd -- silently missing a same-directory
    # binary. Resolve to an absolute path so a relative `binary` still works.
    binary = binary.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    prefix = workdir / "chrM"
    subprocess.run(
        [str(binary), "index", str(fixture_dir / "chrM-human.fa.gz"), str(prefix)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    reads = [str(fixture_dir / "chrM-read_1.fa.gz"), str(fixture_dir / "chrM-read_2.fa.gz")]
    sams: dict[str, Path] = {}
    for label, flags, paired in _MODES:
        out = workdir / f"aln.{label}.sam"
        with out.open("w") as handle:
            subprocess.run(
                [str(binary), "map", *flags, str(prefix), *(reads if paired else reads[:1])],
                stdout=handle,
                stderr=subprocess.DEVNULL,
                check=True,
            )
        sams[label] = out
    return sams


def _modes_are_distinct(stock_digests: dict[str, str]) -> GateResult:
    """Assert every `_MODES` entry actually changes stock's output.

    A mode whose SAM equals another mode's is not a second gate, it is the same
    gate reported twice -- the coverage claim inflates while the evidence does
    not. This fires on the STOCK digests, so it is a property of the mode list and
    upstream's fixture alone; a candidate regression cannot mask it, and cannot
    trigger it either.

    It exists because `-a` was very nearly added as a mode: in `map` it is only
    `flag &= ~MB_F_PAF`, a no-op given SAM is the default, and it was caught by
    running the comparison by hand rather than by anything in this file. The same
    check now runs on every sync, on both architectures.
    """
    by_digest: dict[str, list[str]] = {}
    for label, digest in stock_digests.items():
        by_digest.setdefault(digest, []).append(label)
    duplicates = [labels for labels in by_digest.values() if len(labels) > 1]
    if duplicates:
        groups = "; ".join(" == ".join(labels) for labels in duplicates)
        return GateResult(
            name="modes-are-distinct",
            passed=False,
            detail=f"redundant mode(s) -- identical stock output: {groups}",
        )
    return GateResult(
        name="modes-are-distinct",
        passed=True,
        detail=f"{len(stock_digests)} mode(s), each with distinct stock output",
    )


def run_feature_suites(
    repo: Path, manifest: Manifest, merged: tuple[str, ...] | None
) -> list[GateResult]:
    """Run the `test-*.sh` suites of every merged feature that declares `tests`.

    This is the positive-path coverage `run_gates` deliberately leaves out. Its
    byte-identity comparison proves the merge did not disturb stock behaviour,
    which is the property rerere can silently break -- but it is blind to every
    path a feature turns ON, because those change output on purpose. A feature's
    own suites are the only thing that looks at them, and until this existed
    nothing ran them: `distro-test` covers the Python engine, and the sync
    compiled the aligner without ever invoking `make test`.

    Scoped to `merged` for the same reason the identity coverage line is: a
    dropped feature's scripts are not in the tree, so a bare `make test` step
    would fail with "no rule to make target" rather than a test failure, and
    reporting it as covered would overstate what this build proved.

    Suites run with `repo` as both the working directory and their sole
    argument, matching the `[<minibwa-dir>]` convention the in-tree suites
    already use to locate the binary they exercise.
    """
    results: list[GateResult] = []
    for feature in manifest.features:
        if feature.tests is None:
            continue
        if merged is not None and feature.name not in merged:
            continue
        name = f"tests:{feature.name}"
        directory = repo / feature.tests
        if not directory.is_dir():
            results.append(
                GateResult(
                    name=name,
                    passed=False,
                    detail=f"declared suite directory '{feature.tests}' is not in the build",
                )
            )
            continue
        scripts = sorted(directory.glob("test-*.sh"))
        if not scripts:
            # Zero scripts is a pass over nothing, which is worse than no gate:
            # it reports coverage that does not exist.
            results.append(
                GateResult(name=name, passed=False, detail=f"no test-*.sh in '{feature.tests}'")
            )
            continue
        failures: list[str] = []
        for script in scripts:
            try:
                proc = subprocess.run(
                    ["sh", str(script), str(repo)],
                    cwd=repo,
                    capture_output=True,
                    text=True,
                    timeout=_SUITE_TIMEOUT_S,
                )
                ok, output = proc.returncode == 0, proc.stdout + proc.stderr
            except subprocess.TimeoutExpired:
                ok, output = False, f"timed out after {_SUITE_TIMEOUT_S}s"
            if not ok:
                failures.append(f"{script.name}: {_tail(output)}")
        results.append(
            GateResult(
                name=name,
                passed=not failures,
                detail=(
                    f"{len(scripts)} suite(s) passed"
                    if not failures
                    else f"{len(failures)} of {len(scripts)} failed -- "
                    + " | ".join(failures[:_MAX_REPORTED_FAILURES])
                    + (
                        f" (+{len(failures) - _MAX_REPORTED_FAILURES} more)"
                        if len(failures) > _MAX_REPORTED_FAILURES
                        else ""
                    )
                ),
            )
        )
    return results


def run_gates(
    candidate_bin: Path,
    stock_bin: Path,
    fixture_dir: Path,
    manifest: Manifest,
    workdir: Path,
    merged: tuple[str, ...] | None = None,
    repo: Path | None = None,
) -> list[GateResult]:
    """Run every gate implied by the manifest.

    The comparisons cover all `identical` features and the NEGATIVE case of every
    `conditional` feature at once: with no feature flag set and a reference
    carrying no ALT contigs, none of them may alter a single record -- under any
    of the four upstream invocations in `_MODES`, not merely the default one.

    `merged` names the features the candidate binary actually contains --
    `AssemblyResult.merged`. Optional features drop out routinely, and a gate
    that reads its coverage off the manifest alone credits code that is not in
    the binary it just compared: the one artifact whose job is to be trustworthy
    would overstate what it proved, and disagree with the drop list the same job
    publishes. `None` means "no assembly to scope by, assume the whole manifest".

    Deliberate limitation: the POSITIVE case of a conditional feature (`--meth`,
    `-L`, an ALT index) is not gated here. Those paths change output on purpose,
    so byte-identity says nothing about them; they are covered by each feature
    branch's own tests. This gate proves the merge did not disturb stock behaviour,
    which is the property rerere can silently break.
    """

    def in_build(name: str) -> bool:
        return merged is None or name in merged

    results: list[GateResult] = []
    try:
        stock_sams = align_fixture(stock_bin, fixture_dir, workdir / "stock")
        cand_sams = align_fixture(candidate_bin, fixture_dir, workdir / "cand")
    except (subprocess.CalledProcessError, OSError) as exc:
        results.append(
            GateResult(
                name="default-flags-byte-identity",
                passed=False,
                detail=f"could not run the aligner: {exc}",
            )
        )
    else:
        gated = [f.name for f in manifest.features if f.output in ("identical", "conditional")]
        covered = [name for name in gated if in_build(name)]
        excluded = [name for name in gated if not in_build(name)]
        coverage = f"{len(covered)} feature(s) covered: {', '.join(covered)}"
        if excluded:
            coverage += f"; {len(excluded)} not in this build: {', '.join(excluded)}"
        stock_digests = {label: sam_digest(stock_sams[label]) for label, _, _ in _MODES}
        results.append(_modes_are_distinct(stock_digests))
        for label, flags, paired in _MODES:
            stock = stock_digests[label]
            cand = sam_digest(cand_sams[label])
            passed = stock == cand
            invocation = " ".join(("map", *flags)) + (" (paired)" if paired else " (single-end)")
            results.append(
                GateResult(
                    name=_gate_name(label),
                    passed=passed,
                    # Only the default-flags gate reports coverage; repeating a
                    # fourteen-name list four times buries the one line that
                    # differs when a mode fails.
                    detail=(
                        (coverage if label == "default-flags" else invocation)
                        if passed
                        else f"SAM differs under `{invocation}`: "
                        f"stock {stock[:12]} vs candidate {cand[:12]}"
                    ),
                )
            )

    if repo is not None:
        results.extend(run_feature_suites(repo, manifest, merged))

    for feature in manifest.features:
        # A dropped changes-output feature is not in the binary, so it has
        # nothing to block: failing the release for code this build does not
        # carry is the same misreport in the other direction.
        if feature.output == "changes-output" and in_build(feature.name):
            results.append(
                GateResult(
                    name=f"changes-output:{feature.name}",
                    passed=False,
                    detail="requires a real-data concordance measurement before release",
                )
            )
    return results
