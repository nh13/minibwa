"""Contract tests for the shell the workflows run.

The engine is tested directly; the workflows are not runnable here (no `gh`,
no network, no real git remotes). What can be pinned is the handful of shell
decisions whose failure mode is silent data loss -- exactly where a plausible
one-token edit does the damage.

Every step's script is reached through `yaml_lite.step()`, never through a
whole-file text search: a raw `re.search`/`in` over the file text can be
satisfied by a comment, by a *different* step that happens to share a
substring, or by the step's own name/`if:`/`env:` fields -- none of which run.
Scoping to `step["run"]` rules all three out structurally, and the parser
drops real YAML comments for free (a shell comment *inside* a `run: |` block
is still shell content and is preserved -- see `_executable` below for that
distinction).

For the handful of guards whose whole job is to fail loudly on a bad path,
text presence is not proof: a `case`/`if` that merely still *contains* the
word `exit` can have had that `exit` neutered into an `echo`, and a token
search cannot tell the difference. Those tests instead execute the step's
exact `run:` text with `bash -e -o pipefail` (matching how Actions runs a
`shell: bash` step) against stubbed `git`/`python`/`make`, and assert on the
real exit code and the real argv the stubs observed. That is immune to
restyling (`[ ]` vs `[[ ]]`, quoting, spacing) by construction, and it is
sensitive to exactly the class of mutation that guts a guard while leaving it
looking intact.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

from minibwa_dist.tests import yaml_lite

WORKFLOWS = Path(__file__).resolve().parents[2] / ".github" / "workflows"


def _parsed(name: str) -> dict[str, object]:
    return yaml_lite.parse((WORKFLOWS / name).read_text())


def _sync() -> dict[str, object]:
    return _parsed("distro-sync.yml")


def _release() -> dict[str, object]:
    return _parsed("distro-release.yml")


def _reconcile() -> dict[str, object]:
    return _parsed("distro-reconcile.yml")


def _ship() -> dict[str, object]:
    return _parsed("distro-ship.yml")


def _test_workflow() -> dict[str, object]:
    return _parsed("distro-test.yml")


def _run(workflow: dict[str, object], job: str, name: str) -> str:
    """The exact `run:` script text of the step named `name` in `job`."""
    value = yaml_lite.step(workflow, job, name=name).get("run")
    assert isinstance(value, str), f"step {name!r} in job {job!r} has no run: script"
    return value


def _all_run_text(workflow: dict[str, object]) -> str:
    """Every step's `run:` text in `workflow`, across every job, concatenated.

    For assertions that are legitimately whole-workflow ("this pattern must
    appear nowhere"), not whole-*file* -- `name:`/`if:`/comment text is still
    excluded, only `run:` bodies are joined.
    """
    jobs = workflow.get("jobs")
    assert isinstance(jobs, dict)
    chunks: list[str] = []
    for job in jobs.values():
        for entry in job.get("steps", []):
            run = entry.get("run") if isinstance(entry, dict) else None
            if isinstance(run, str):
                chunks.append(run)
    return "\n".join(chunks)


def _mktemp_path() -> Path:
    """A fresh, empty file path a stub can log to, closed but not deleted.

    `tempfile.mkstemp` returns an open file descriptor as well as the path;
    the descriptor is only needed to guarantee the name is unique and unused,
    so it is closed immediately -- the stub subprocess writes to the path via
    a fresh `>>` open, not through this handle.
    """
    fd, path = tempfile.mkstemp()
    os.close(fd)
    return Path(path)


def _executable(text: str) -> str:
    """`text` with shell comment lines dropped, for "this command is gone"
    assertions against a step's own `run:` body.

    A `#`-led line here is a genuine shell comment (the YAML layer already
    stripped anything that was a YAML comment when `yaml_lite` parsed the
    block scalar) -- so this is "what bash actually executes", not "what
    appears in the script's text". That distinction is what defeats a
    mutation that reintroduces a removed command as a comment sitting next to
    an innocuous replacement.
    """
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))


# -- behavioural harness: run a step's script for real, against stub tools --


def _sh(
    script: str,
    *,
    stubs: dict[str, str],
    cwd_executables: dict[str, str] | None = None,
    cwd_files: dict[str, str] | None = None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run `script` -- one step's exact `run:` text -- the way Actions' own
    `shell: bash` executes it, with `stubs` (command name -> script body)
    standing in for external tools on `PATH`.

    `cwd_executables` covers `./minibwa`-style calls, which resolve relative
    to the working directory and are never found via `PATH`. `PATH` is the
    inherited one with the stub directory prepended, not a hardcoded list --
    this suite runs on both a developer's machine and `ubuntu-latest`, and a
    hardcoded path list would quietly stop finding `bash`/`git` on one of them.
    """
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        bin_dir = work / "bin"
        bin_dir.mkdir()
        for command, body in stubs.items():
            stub = bin_dir / command
            stub.write_text(f"#!/bin/bash\n{body}\n")
            stub.chmod(0o755)
        for name, body in (cwd_executables or {}).items():
            exe = work / name
            exe.write_text(f"#!/bin/bash\n{body}\n")
            exe.chmod(0o755)
        for name, content in (cwd_files or {}).items():
            (work / name).write_text(content)
        env = {
            **os.environ,
            "PATH": f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}",
            **(extra_env or {}),
        }
        return subprocess.run(
            ["bash", "-e", "-o", "pipefail", "-c", script],
            cwd=work,
            env=env,
            capture_output=True,
            text=True,
        )


def test_the_rerere_cache_fetch_distinguishes_absent_from_failed() -> None:
    """`|| true` on this fetch made a network or auth failure look like a first
    run, and the save step answers a first run by building a fresh orphan.
    Dropping `|| true` outright is not the fix either -- fetching a genuinely
    nonexistent ref exits 128, which would break every cold start. The remote
    has to be asked first, by exit code: exit 2 (no such ref) must fall
    through to a normal run, anything else must take the job down.

    Mutation this guards: the `[ "$rc" -ne 2 ]` branch stops exiting and just
    echoes -- a real fetch failure (rc=128, say) would then be swallowed and
    the job would carry on as if this were a first run.
    """
    script = _run(_sync(), "assemble", "Fetch the fork, then upstream")
    git_stub = """
    case "$1" in
      fetch) exit 0 ;;
      ls-remote) exit "${LS_REMOTE_EXIT:-0}" ;;
      remote) exit 0 ;;
      push) exit 0 ;;
      *) echo "unexpected git $*" >&2; exit 99 ;;
    esac
    """

    failure = _sh(script, stubs={"git": git_stub}, extra_env={"LS_REMOTE_EXIT": "128"})
    assert failure.returncode == 128, (
        f"a real fetch failure (exit 128) must take the job down, got {failure.returncode}: "
        f"{failure.stderr}"
    )

    absent = _sh(script, stubs={"git": git_stub}, extra_env={"LS_REMOTE_EXIT": "2"})
    assert absent.returncode == 0, (
        f"exit 2 (no such ref) is a genuine first run and must not fail the job: {absent.stderr}"
    )


def test_the_rerere_cache_is_never_force_pushed() -> None:
    """A fresh orphan shares no history with the published branch, so a plain
    push is rejected as non-fast-forward -- git itself refuses the wipe.
    `--force` turns that refusal into the wipe, and `|| true` hides the
    rejection (or any other push failure) from the job.

    Mutation this guards: `git push origin rerere-cache` becomes
    `git push --force origin rerere-cache || true`, with an innocuous comment
    left above mentioning the plain form -- a comment that a whole-file
    substring search over "git push...rerere-cache" can match instead of the
    real, mutated line. Executing the real line sidesteps that: the comment
    never runs, so it cannot appear in the stub's logged argv, and a `|| true`
    is caught directly by the exit code rather than by finding the token.
    """
    script = _run(_sync(), "assemble", "Save the rerere cache")
    log = _mktemp_path()
    git_stub = f"""
    case "$1" in
      rev-parse) exit 0 ;;
      push) echo "$*" >> "{log}"; exit 17 ;;
      *) echo "unexpected git $*" >&2; exit 99 ;;
    esac
    """

    result = _sh(script, stubs={"git": git_stub})

    assert result.returncode == 17, (
        "a real push failure must fail the step, not be swallowed by `|| true`: "
        f"rc={result.returncode}, stderr={result.stderr}"
    )
    pushed = log.read_text() if log.exists() else ""
    log.unlink(missing_ok=True)
    assert "push" in pushed, "the guard's `then` branch must have pushed at all"
    assert "--force" not in pushed, f"the cache branch must never be force-pushed: {pushed!r}"


def test_dist_next_is_still_force_pushed() -> None:
    """The non-force rule is about the cache only: `dist-next` is reassembled
    from scratch every run and is meant to be replaced.

    Token-based, not a literal substring: flag order (`--force origin` vs.
    `origin --force`) is not the thing under test, `--force` actually being
    there is.
    """
    text = _run(_sync(), "assemble", "Push dist-next and open the review PR")
    line = next(line for line in text.splitlines() if "dist-next" in line and "git push" in line)
    tokens = line.split()
    assert "--force" in tokens and "dist-next" in tokens, line


def test_the_release_lineage_guard_compares_the_recorded_assembly_base() -> None:
    """`git merge-base --is-ancestor "$TAG" HEAD` used to stand here, and it
    passes for exactly the case the step's own comment says it prevents:
    `dist` is reassembled nightly from upstream master tip, so any dist built
    after v0.6 still contains v0.6. Ancestry cannot express "assembled at this
    tag"; an equal sha can.

    Mutation this guards: the mismatch branch's `exit 1` is dropped, so a
    `dist` assembled from the wrong base would release anyway.
    """
    step = yaml_lite.step(
        _release(), "release", name="Confirm dist was assembled from this exact upstream tag"
    )
    text = step.get("run")
    assert isinstance(text, str)
    assert "merge-base --is-ancestor" not in _executable(text)
    assert "jq -r '.base // empty' dist-manifest.json" in text
    assert 'git rev-parse "${TAG}^{commit}"' in text
    assert '[ "$BASE" = "$TAG_SHA" ]' in text

    git_stub = """
    case "$1" in
      fetch) exit 0 ;;
      remote) exit 0 ;;
      rev-parse) echo "${TAG_SHA_STUB:?}"; exit 0 ;;
      *) echo "unexpected git $*" >&2; exit 99 ;;
    esac
    """
    result = _sh(
        text,
        stubs={"git": git_stub},
        cwd_files={"dist-manifest.json": json.dumps({"base": "aaaaaaaaaaaa"})},
        extra_env={"TAG": "v0.6", "TAG_SHA_STUB": "bbbbbbbbbbbb"},
    )
    assert result.returncode == 1, (
        f"a base/tag mismatch must fail the release: rc={result.returncode}, {result.stderr}"
    )


def test_the_release_gates_the_tooling_the_tag_actually_carries() -> None:
    """Overlaying origin/ops/distro's tip made the gate and the release notes
    describe a manifest the tag does not contain -- and the `git checkout --
    minibwa_dist pyproject.toml` meant to undo it restored from the index that same
    overlay had already written, so it undid nothing.
    """
    text = _executable(_all_run_text(_release()))
    assert "git checkout origin/ops/distro" not in text
    assert "git checkout -- minibwa_dist pyproject.toml" not in text


def test_the_sync_gate_step_scopes_coverage_and_propagates_a_failure() -> None:
    """`--assembly` scopes the reported coverage to what actually merged: an
    optional feature that conflicted is not in this binary, and a gate that
    names it as covered contradicts the drop list this same job publishes.
    And the aggregation that turns a gate failure into a nonzero exit is the
    whole point of the step -- it used to live as a hand-copied heredoc,
    reachable by neither `ruff check minibwa_dist/` nor `pytest minibwa_dist/tests`.

    Mutations this guards, both on the "Output gates" step: `--assembly`
    moved out of the actual command into a step comment (so the real
    invocation stops scoping coverage); and `python -m minibwa_dist.cli gates`
    gaining a trailing `|| true` (so a real gate failure stops failing the
    job). A logged-argv check catches the first even though the flag is still
    *somewhere* in the file; an exit-code check catches the second even
    though the command is still spelled out in full.
    """
    script = _run(_sync(), "assemble", "Output gates")
    assert "python -m minibwa_dist.cli gates" in script
    assert "from minibwa_dist.gates import run_gates" not in _executable(script)

    log = _mktemp_path()
    stubs = {
        "make": "exit 0",
        "git": 'case "$1" in worktree) exit 0 ;; *) echo "unexpected git $*" >&2; exit 99 ;; esac',
        "python": f"""
        echo "$*" >> "{log}"
        case "$*" in
          *"minibwa_dist.cli gates"*) exit 1 ;;
          *) exit 0 ;;
        esac
        """,
    }

    result = _sh(script, stubs=stubs)

    assert result.returncode == 1, (
        f"a gate failure must fail the step, not be swallowed: rc={result.returncode}"
    )
    invocation = next(
        (line for line in log.read_text().splitlines() if "minibwa_dist.cli gates" in line), ""
    )
    log.unlink(missing_ok=True)
    assert "--assembly /tmp/assembly.json" in invocation, (
        f"the real invocation must carry --assembly, not just the file text: {invocation!r}"
    )


def test_the_release_gate_step_scopes_coverage_and_propagates_a_failure() -> None:
    """The release job's own copy of the same gate step -- see
    `test_the_sync_gate_step_scopes_coverage_and_propagates_a_failure` for why
    both the argument and the exit code are checked by execution rather than
    by text search.
    """
    script = _run(_release(), "release", "Build, gate and verify BEFORE tagging")
    assert "python -m minibwa_dist.cli gates" in script
    assert "from minibwa_dist.gates import run_gates" not in _executable(script)

    log = _mktemp_path()
    stubs = {
        "make": "exit 0",
        "git": 'case "$1" in worktree) exit 0 ;; *) echo "unexpected git $*" >&2; exit 99 ;; esac',
        "python": f"""
        echo "$*" >> "{log}"
        case "$*" in
          *"minibwa_dist.cli gates"*) exit 1 ;;
          *) exit 0 ;;
        esac
        """,
    }

    result = _sh(
        script,
        stubs=stubs,
        cwd_executables={"minibwa": 'echo "$VERSION"'},
        extra_env={"VERSION": "9.9.9-test", "TAG": "v0.6"},
    )

    assert result.returncode == 1, (
        f"a gate failure must fail the step, not be swallowed: rc={result.returncode}"
    )
    invocation = next(
        (line for line in log.read_text().splitlines() if "minibwa_dist.cli gates" in line), ""
    )
    log.unlink(missing_ok=True)
    assert "--assembly dist-manifest.json" in invocation, (
        f"the real invocation must carry --assembly, not just the file text: {invocation!r}"
    )


def test_the_release_stamps_the_version_through_the_cli_not_sed() -> None:
    """Two implementations of the same header rewrite, with different match
    semantics and only one of them tested -- and the sed one spliced an
    operator-supplied version into a substitution as syntax.
    """
    text = _run(_release(), "release", "Stamp the version")
    assert "python -m minibwa_dist.cli stamp" in text
    assert "sed -i" not in _executable(text)


def test_the_release_validates_its_dispatch_inputs() -> None:
    """workflow_dispatch inputs are attacker-controllable text; the job hoists
    them to env to keep them out of the Actions template, then feeds them to a
    version string, a filename and `git tag`.

    Mutation this guards: both `case` arms stop exiting and just echo instead
    -- so a `$REV`/`$TAG` that fails validation would sail through to `git
    tag` anyway. Each arm is exercised on its own bad input, with the other
    input held valid, so a mutation in either one is caught on its own.
    """
    script = _run(_release(), "release", "Validate the dispatch inputs")
    # Whitespace-tolerant, not a literal substring: incidental spacing around
    # `in` is not the thing under test, the patterns themselves are.
    assert re.search(r"case \"\$REV\" in\s+''\|\*\[!0-9\]\*\)", script)
    assert re.search(r"case \"\$TAG\" in\s+''\|\*\[!A-Za-z0-9\._-\]\*\)", script)

    bad_rev = _sh(script, stubs={}, extra_env={"REV": "not-a-number", "TAG": "v0.6"})
    assert bad_rev.returncode != 0, "a non-numeric revision must fail validation"

    bad_tag = _sh(script, stubs={}, extra_env={"REV": "1", "TAG": "v0.6; rm -rf /"})
    assert bad_tag.returncode != 0, "a tag with disallowed characters must fail validation"

    both_valid = _sh(script, stubs={}, extra_env={"REV": "1", "TAG": "v0.6"})
    assert both_valid.returncode == 0, "valid inputs must not be rejected"


def test_the_nominations_issue_title_carries_no_varying_count() -> None:
    """The dedupe searches by title, so a count in it means a changed count
    files a duplicate on top of the still-open issue -- the failure mode the
    sibling step in distro-sync.yml explicitly reasons its way out of.
    """
    text = _run(_reconcile(), "reconcile", "Monthly upstreaming nominations")
    assert 'TITLE="upstream: features never offered to lh3/minibwa"' in text
    assert "$N feature(s)" not in _executable(text)
    assert "gh issue edit" in text, "a dedupe hit must refresh the body, not drop the update"


def test_the_reconciler_files_an_issue_when_it_fails() -> None:
    """Without this the reconciler goes quiet rather than loud, and stale
    statuses are noticed only by someone reading the manifest.
    """
    step = yaml_lite.step(_reconcile(), "reconcile", name="File an issue on failure")
    assert step.get("if") == "failure()"
    text = step.get("run")
    assert isinstance(text, str)
    assert 'TITLE="reconcile: run failed"' in text


def test_the_workflows_do_not_hand_parse_the_manifest() -> None:
    """Schema knowledge belongs in `minibwa_dist.manifest`. Two heredocs used to dig
    through the raw TOML, and had already drifted apart on what a manifest
    with no [[feature]] table even means.
    """
    for workflow in (_sync(), _reconcile()):
        assert "tomllib" not in _executable(_all_run_text(workflow))
    assert "from minibwa_dist.manifest import load_manifest" in _run(
        _reconcile(), "reconcile", "Collect live PR states"
    )


def test_the_nochange_guard_ignores_only_derived_artifacts() -> None:
    """A no-op sync must not open a PR.

    MB_VERSION embeds the assembly sha and dist-manifest.json records it, so both
    differ on every rebuild. If the "nothing to ship" check compared raw trees it
    could never fire and the bot would file a PR every day, which is how a review
    channel stops being read. Assert the guard excludes exactly those two and
    nothing else -- excluding minibwa.h wholesale would hide real feature changes.
    """
    script = _run(_sync(), "assemble", "Re-parent dist-next onto dist")
    assert "dist-manifest.json" in script
    assert "#define MB_VERSION" in script
    assert "':!minibwa.h'" not in script, "excluding minibwa.h wholesale would hide feature changes"
    assert "nochange=1" in script


def test_the_ship_step_survives_github_marking_the_pr_merged_itself() -> None:
    """Force-updating `dist` to this PR's head makes that head an ancestor of
    the base branch, so GitHub marks the PR MERGED on its own -- asynchronously,
    so whether that lands before or after this step is a race. `gh pr close`
    errors on an already-merged PR, so the ship run went red *after* the push
    that was its entire job had succeeded. A red run on a good ship is how
    people learn to stop reading ship runs.

    The postcondition asserted here is "this PR is not left open", never "this
    step called `gh pr close`" -- checking the call would reject a perfectly
    good attempt-then-verify rewrite, and the end state is the thing that
    matters. The third arm is what carries the test: the obvious fix, a bare
    `gh pr close || true`, passes the first two and fails only that one,
    because it would swallow a genuine auth or API failure and leave the PR
    open with the run still green.
    """
    script = _run(_ship(), "ship", "Close the PR")
    state_file = _mktemp_path()
    log = _mktemp_path()
    gh_stub = f"""
    echo "$*" >> "{log}"
    case "$1 $2" in
      "pr comment") exit 0 ;;
      "pr view") cat "{state_file}" ;;
      "pr close")
        rc="${{CLOSE_EXIT:-0}}"
        [ "$rc" -eq 0 ] && echo CLOSED > "{state_file}"
        exit "$rc" ;;
      *) echo "unexpected gh $*" >&2; exit 99 ;;
    esac
    """

    def ship(initial_state: str, close_exit: str = "0") -> tuple[int, str]:
        state_file.write_text(f"{initial_state}\n")
        log.write_text("")
        result = _sh(
            script,
            stubs={"gh": gh_stub},
            extra_env={
                "PR_NUMBER": "17",
                "HEAD_SHA": "800583d3fa9a6658c4a4f3254070929886163ac8",
                "CLOSE_EXIT": close_exit,
            },
        )
        return result.returncode, log.read_text()

    already_merged_rc, _ = ship("MERGED")
    assert already_merged_rc == 0, (
        "GitHub auto-merging the PR is the expected outcome of the force-push, not a failure: "
        f"rc={already_merged_rc}"
    )

    open_rc, open_log = ship("OPEN")
    assert open_rc == 0, f"a still-open PR must be closed cleanly: rc={open_rc}"
    assert "pr close" in open_log, (
        f"the race can go the other way -- an open PR must still be closed: {open_log!r}"
    )

    stuck_rc, _ = ship("OPEN", close_exit="1")
    assert stuck_rc != 0, (
        "a close that fails and leaves the PR open must fail the step -- this is the arm a bare "
        f"`|| true` gets wrong: rc={stuck_rc}"
    )

    state_file.unlink(missing_ok=True)
    log.unlink(missing_ok=True)


def test_the_test_workflow_pins_its_python_tooling() -> None:
    """Every `uses:` in these workflows is pinned to a SHA; an unpinned linter
    is the same class of "CI turns red for something that is not the change".
    """
    text = _all_run_text(_test_workflow())
    assert re.search(r"pipx install ruff==\d", text)
    assert re.search(r"pip install pytest==\d", text)


def test_the_dropped_feature_issue_refreshes_on_a_dedupe_hit() -> None:
    """A stable title means the dedupe hits whenever *any* drop issue is open.

    The step used to `exit 0` on a hit, so a drop of a different feature set was
    swallowed by an issue naming an unrelated one -- which is how a real drop
    reached nobody. The reconciler already edits on a hit; this step must too,
    and it must refresh the title as well, because the title names the set.
    """
    text = _run(_sync(), "assemble", "File an issue for dropped features")
    assert "gh issue edit" in text, "a dedupe hit must refresh the issue, not drop the update"
    assert "exit 0" not in _executable(text), "a dedupe hit must not silently skip the update"
    assert "--title" in text, "the title names the dropped set, so a refresh must update it"


def test_the_sync_passes_a_regression_baseline() -> None:
    """Without it, a feature that merged yesterday and does not today exits zero.

    `origin/dist-next` still holds the previous build at the moment sync runs --
    this run force-pushes it only afterwards -- so it is the ref that answers
    "what did the last build contain".
    """
    text = _run(_sync(), "assemble", "Assemble dist-next")
    assert "--regression-baseline origin/dist-next" in _executable(text)


def test_the_assembly_json_reaches_the_log_even_when_the_sync_fails() -> None:
    """A regression exits nonzero, and `bash -e` would abort before the `cat`.

    The run would then go red with the reason only on stderr and the assembly
    itself -- which features dropped, and what they conflict in -- never shown.
    Capture the status, print, then exit with it, so a red run is diagnosable
    from the log alone. The exit code must still propagate: masking it is what
    the step's own "redirect, do not pipe" note already guards against.
    """
    text = _executable(_run(_sync(), "assemble", "Assemble dist-next"))
    assert "|| rc=$?" in text, "the sync's exit status must be captured, not aborted on"
    assert "exit $rc" in text, "the captured status must still fail the step"
    assert text.index("cat /tmp/assembly.json") < text.index("exit $rc")


def test_a_hard_assembly_failure_exits_before_the_jq_parsing() -> None:
    """A required-feature conflict prints no JSON, so the file is empty.

    Falling through to `jq` would then fail the step on a parse error, burying
    the real cause -- which is on stderr -- under a confusing one. A regression
    is different: it DOES emit valid JSON, so it must still reach the outputs.
    """
    text = _executable(_run(_sync(), "assemble", "Assemble dist-next"))
    guard = "jq -e . /tmp/assembly.json"
    assert guard in text, "a nonzero run with no parsable JSON must exit before jq"
    assert text.index(guard) < text.index("dropped=$(jq")


def test_the_gates_step_builds_the_api_probes_it_depends_on() -> None:
    """Feature suites shell out to the api-test probes, so the gates step must
    build them rather than inherit them from an earlier step's side effect.
    Without this the gate fails with "probe not built" -- a build error dressed
    up as a test failure -- if the steps are ever reordered or split.
    """
    text = _executable(_run(_sync(), "assemble", "Output gates"))
    assert "make -C api-test" in text or "-C api-test" in text


def test_both_gate_jobs_point_repo_at_the_assembled_tree() -> None:
    """`--repo` selects where the feature suites are run, and the two jobs lay
    their workspaces out differently: the x86 job runs from the assembled
    checkout, while the arm64 job checks out `ops/distro` and builds the
    candidate in a worktree. A default would silently aim the arm64 gate at the
    tooling checkout, which carries no suites at all.
    """
    x86 = _executable(_run(_sync(), "assemble", "Output gates"))
    assert "--repo ." in x86

    arm = _executable(_run(_sync(), "gate-arm64", "Re-gate the published assembly on arm64"))
    assert "--repo /tmp/cand" in arm, "the arm64 gate must run the suites in the candidate tree"
    assert "/tmp/cand/api-test" in arm, "and must build the probes those suites shell out to"


def test_the_drop_issue_edit_targets_an_exactly_shaped_title() -> None:
    """`gh issue edit` rewrites title AND body, so the target cannot be chosen by
    a token search alone -- `sync: dropped in:title` matches on two ordinary
    words, and `.[0]` would hand the bot whatever came back first. The previous
    code was safe only because a false hit merely skipped the update.
    """
    text = _executable(_run(_sync(), "assemble", "File an issue for dropped features"))
    assert 'startswith("sync: dropped ")' in text
    assert 'endswith(" -- would not merge")' in text


def test_the_drop_issue_step_runs_even_when_the_assembly_failed() -> None:
    """A regression exits the assemble step nonzero, and the implicit `success()`
    would skip the one issue that names which features dropped. That step writes
    its outputs before exiting precisely so this one can still consume them.
    """
    step = yaml_lite.step(_sync(), "assemble", name="File an issue for dropped features")
    condition = str(step.get("if"))
    assert condition.startswith("always()"), f"expected always() guard, got {condition!r}"
    assert "dropped != '[]'" in condition, "still only when something actually dropped"


def test_the_suite_timeout_fits_inside_the_job_timeout() -> None:
    """A per-suite budget the whole set can outlast is a guard that never fires:
    the runner kills the job first and prints nothing.
    """
    from minibwa_dist.gates import _SUITE_TIMEOUT_S

    jobs = _sync()["jobs"]
    caps = [int(job["timeout-minutes"]) for job in jobs.values() if "timeout-minutes" in job]
    assert caps, "the gate jobs must declare timeout-minutes for this bound to mean anything"
    # 15 in-tree ALT suites today; leave room for a set that grows.
    assert _SUITE_TIMEOUT_S * 20 <= min(caps) * 60, (
        f"{_SUITE_TIMEOUT_S}s x 20 suites exceeds the {min(caps)}min job cap"
    )
