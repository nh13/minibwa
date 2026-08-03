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


def align_fixture(binary: Path, fixture_dir: Path, workdir: Path) -> Path:
    """Index the chrM fixture with `binary` and align the read pair with it.

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
    out = workdir / "aln.sam"
    with out.open("w") as handle:
        subprocess.run(
            [
                str(binary),
                "map",
                str(prefix),
                str(fixture_dir / "chrM-read_1.fa.gz"),
                str(fixture_dir / "chrM-read_2.fa.gz"),
            ],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    return out


def run_gates(
    candidate_bin: Path,
    stock_bin: Path,
    fixture_dir: Path,
    manifest: Manifest,
    workdir: Path,
    merged: tuple[str, ...] | None = None,
) -> list[GateResult]:
    """Run every gate implied by the manifest.

    One default-flags comparison covers all `identical` features and the NEGATIVE
    case of every `conditional` feature at once: with no feature flag set and a
    reference carrying no ALT contigs, none of them may alter a single record.

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
        stock = sam_digest(align_fixture(stock_bin, fixture_dir, workdir / "stock"))
        cand = sam_digest(align_fixture(candidate_bin, fixture_dir, workdir / "cand"))
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
        detail = f"{len(covered)} feature(s) covered: {', '.join(covered)}"
        if excluded:
            detail += f"; {len(excluded)} not in this build: {', '.join(excluded)}"
        passed = stock == cand
        results.append(
            GateResult(
                name="default-flags-byte-identity",
                passed=passed,
                detail=(
                    detail
                    if passed
                    else f"SAM differs: stock {stock[:12]} vs candidate {cand[:12]}"
                ),
            )
        )

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
