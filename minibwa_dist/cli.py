"""Command-line entry points for the distribution workflows.

Kept thin on purpose: every non-trivial decision lives in a tested module.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from minibwa_dist import rerere_cache
from minibwa_dist.assemble import (
    AssemblyError,
    TrainingStop,
    assemble,
    previously_merged,
    regressions,
    set_version,
    stamp_dev_version,
)
from minibwa_dist.gates import run_gates
from minibwa_dist.gitutil import GitError
from minibwa_dist.manifest import ManifestError, load_manifest, parse_manifest
from minibwa_dist.reconcile import apply_changes, graduated, nomination_candidates, reconcile
from minibwa_dist.render import render_graveyard, render_release_notes, update_readme
from minibwa_dist.rerere_cache import OrphanOverPublishedCache


def _sync(args: argparse.Namespace) -> int:
    repo = Path(args.repo)
    manifest = load_manifest(Path(args.manifest))
    restored = rerere_cache.restore(repo)
    print(f"rerere: restored {restored} resolution(s)", file=sys.stderr)
    try:
        result = assemble(
            repo, manifest, args.base, args.out, remote=args.remote, stop_at=args.stop_at
        )
    except TrainingStop:
        # `--stop-at` left a conflict in the worktree on purpose (see
        # assemble()'s docstring). That is the expected, successful outcome of
        # training mode, not a build failure, so it must not print "FAILED" or
        # exit nonzero. A genuine required-feature conflict can still fire while
        # `--stop-at` is set -- it raises the base class, lands in the handler
        # below, and still exits 1.
        #
        # The resolution is only worth recording if it reaches the fork: `save()`
        # ran in the `finally` below at the moment of the conflict, so it has the
        # pre-image and nothing else. Spell out the whole follow-up, naming the
        # remote the caller already chose -- `origin` is upstream in a local
        # checkout, and a resolution pushed there goes to the wrong repository.
        print(
            f"stopped at '{args.stop_at}' for rerere training -- resolve, commit, "
            "then `git reset --hard HEAD~1`. Re-run this same sync to record the "
            f"resolution on the local `{rerere_cache.CACHE_BRANCH}` branch, then "
            f"`git push {args.remote} {rerere_cache.CACHE_BRANCH}` to publish it.",
            file=sys.stderr,
        )
        return 0
    except AssemblyError as exc:
        print(f"ASSEMBLY FAILED: {exc}", file=sys.stderr)
        return 1
    finally:
        # Save even on failure: a partial run may have recorded new pre-images,
        # and throwing them away means paying for the same conflicts next time.
        rerere_cache.save(repo)

    version = stamp_dev_version(repo, args.out) if args.stamp else None
    regressed = (
        regressions(previously_merged(repo, args.regression_baseline), manifest, result.merged)
        if args.regression_baseline
        else ()
    )
    # Emit the JSON even when regressed: the caller redirects this to a file and
    # reads it after the step fails, and a run that fails silently is worse than
    # the drop it is reporting.
    print(
        json.dumps(
            {
                "base": result.base,
                "head": result.head,
                "version": version,
                "merged": list(result.merged),
                "dropped": [
                    {"name": d.name, "files": list(d.conflicting_files)} for d in result.dropped
                ],
                "regressed": list(regressed),
            }
        )
    )
    if regressed:
        print(
            f"REGRESSION: {', '.join(regressed)} merged in the build at "
            f"'{args.regression_baseline}' and no longer merge. Resolve the conflict with "
            f"`--stop-at <name>` and publish the resolution to `{rerere_cache.CACHE_BRANCH}`, "
            "or retire the feature from the manifest if the drop is intended.",
            file=sys.stderr,
        )
        return 1
    return 0


def _gates(args: argparse.Namespace) -> int:
    """Run the output gates and turn their verdicts into an exit code.

    A subcommand rather than two hand-copied heredocs in two workflows: the
    aggregation (`all(...)`) is the whole point of the step, and neither copy
    was reachable by `ruff check minibwa_dist/` or `pytest minibwa_dist/tests`.

    `--repo` defaults to the working directory because the gates step already
    runs from the assembled checkout; the feature suites need to find the binary
    they were just built alongside.
    """
    manifest = load_manifest(Path(args.manifest))
    merged: tuple[str, ...] | None = None
    if args.assembly:
        merged = tuple(json.loads(Path(args.assembly).read_text())["merged"])
    results = run_gates(
        Path(args.candidate),
        Path(args.stock),
        Path(args.fixtures),
        manifest,
        Path(args.workdir),
        merged=merged,
        repo=Path(args.repo),
    )
    for result in results:
        print(f"{'PASS' if result.passed else 'FAIL'}  {result.name}: {result.detail}")
    return 0 if all(result.passed for result in results) else 1


def _stamp(args: argparse.Namespace) -> int:
    """Set MB_VERSION to an exact version string, for the release workflow."""
    previous = set_version(Path(args.repo), args.version)
    print(f"MB_VERSION: {previous} -> {args.version}", file=sys.stderr)
    return 0


def _reconcile(args: argparse.Namespace) -> int:
    manifest_path = Path(args.manifest)
    manifest = load_manifest(manifest_path)
    pr_states = json.loads(Path(args.pr_states).read_text())
    changes = reconcile(manifest, {int(k): v for k, v in pr_states.items()})
    if changes and args.write:
        new_text = apply_changes(manifest_path.read_text(), changes)
        # Refuse to write a manifest this same engine cannot load back: CI's
        # own load of the committed file is not the first opportunity to
        # notice a malformed rewrite, this is.
        parse_manifest(new_text)
        manifest_path.write_text(new_text)
    print(
        json.dumps(
            {
                "changes": [{"feature": c.feature, "old": c.old, "new": c.new} for c in changes],
                "graduated": graduated(manifest, changes),
                "nominations": nomination_candidates(manifest),
            }
        )
    )
    return 0


def _notes(args: argparse.Namespace) -> int:
    manifest = load_manifest(Path(args.manifest))
    dropped = tuple(args.dropped) if args.dropped else ()
    print(
        render_release_notes(
            manifest, Path(args.upstream_notes).read_text(), args.upstream_tag, dropped
        )
    )
    return 0


def _docs(args: argparse.Namespace) -> int:
    """Regenerate README.md's distribution block and GRAVEYARD.md in place."""
    manifest = load_manifest(Path(args.manifest))
    update_readme(Path(args.readme), manifest)
    Path(args.graveyard).write_text(render_graveyard(manifest) + "\n")
    print(f"rendered {args.readme} and {args.graveyard}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="minibwa_dist")
    parser.add_argument(
        "--manifest", default="minibwa_dist/features.toml", help="path to the distribution manifest"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sync", help="assemble dist-next")
    p.add_argument("--repo", default=".", help="path to the minibwa checkout to assemble in")
    p.add_argument(
        "--base", default="master", help="ref to assemble onto, normally upstream master"
    )
    p.add_argument("--out", default="dist-next", help="branch to (re)create with the assembly")
    p.add_argument(
        "--remote",
        required=True,
        help="remote namespace for feature branches: nh13 locally, origin in CI",
    )
    p.add_argument(
        "--stamp",
        action="store_true",
        help="rewrite MB_VERSION so the build identifies itself as downstream",
    )
    p.add_argument(
        "--regression-baseline",
        default=None,
        metavar="REF",
        help="ref whose dist-manifest.json records the previous build; a feature that "
        "merged there and does not now fails the run",
    )
    p.add_argument(
        "--stop-at",
        default=None,
        help="training mode: leave this feature's conflict in the worktree",
    )
    p.set_defaults(func=_sync)

    p = sub.add_parser("reconcile", help="compare recorded statuses to live PR states")
    p.add_argument("--pr-states", required=True, help="JSON file: {pr: state}")
    p.add_argument("--write", action="store_true", help="rewrite features.toml in place")
    p.set_defaults(func=_reconcile)

    p = sub.add_parser("notes", help="render release notes")
    p.add_argument("--upstream-notes", required=True, help="file holding upstream's own notes")
    p.add_argument("--upstream-tag", required=True, help="upstream tag being mirrored, e.g. v0.6")
    p.add_argument(
        "--dropped", nargs="*", help="feature names excluded from this build, space separated"
    )
    p.set_defaults(func=_notes)

    p = sub.add_parser("docs", help="regenerate README.md block and GRAVEYARD.md")
    p.add_argument("--readme", default="README.md", help="README to rewrite the distro block in")
    p.add_argument("--graveyard", default="GRAVEYARD.md", help="graveyard file to overwrite")
    p.set_defaults(func=_docs)

    p = sub.add_parser("gates", help="run the output gates; exit 1 if any fails")
    p.add_argument("--candidate", required=True, help="the distribution binary under test")
    p.add_argument("--stock", required=True, help="a binary built from unmodified upstream")
    p.add_argument("--fixtures", required=True, help="directory holding upstream's chrM fixtures")
    p.add_argument("--workdir", required=True, help="scratch directory for indexes and SAMs")
    p.add_argument(
        "--assembly",
        default=None,
        help="assembly JSON, so coverage names only the features this binary contains",
    )
    p.add_argument(
        "--repo",
        default=".",
        help="the assembled checkout; each merged feature's declared test suites are run in it",
    )
    p.set_defaults(func=_gates)

    p = sub.add_parser("stamp", help="set MB_VERSION to an exact version string")
    p.add_argument("--repo", default=".", help="path to the minibwa checkout to stamp")
    p.add_argument("--version", required=True, help="version to write, without a leading 'v'")
    p.set_defaults(func=_stamp)

    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ManifestError, GitError, OrphanOverPublishedCache) as exc:
        # A malformed manifest, a failed git call, and a rerere cache that is
        # published but was never fetched locally are all ordinary operating
        # conditions for an unattended workflow -- the last of those is what a
        # fresh `git clone` looks like. Report them as one line naming the
        # cause, not as a traceback that discards an assembly that actually
        # succeeded (`_sync` saves the cache from its `finally`, so this can
        # fire after `sync` has already done its real work).
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
