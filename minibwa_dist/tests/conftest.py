"""Synthetic git repositories for testing the assembly engine.

Test data is generated, never committed — see the project conventions.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def run(repo: Path, *args: str) -> str:
    """Run a git command in `repo` and return stdout."""
    result = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def commit_file(repo: Path, name: str, content: str, message: str) -> str:
    """Write `name`, commit it, and return the new commit sha."""
    (repo / name).write_text(content)
    run(repo, "add", name)
    run(repo, "-c", "commit.gpgsign=false", "commit", "-m", message)
    return run(repo, "rev-parse", "HEAD")


def merge_and_resolve(repo: Path, branch: str, path: str, resolution: str) -> None:
    """Merge `branch`, resolve the conflict in `path` to `resolution`, and commit.

    The merge is expected to conflict, so its nonzero exit is not an error —
    that is the whole point. Enables `rerere.enabled` itself before merging:
    git only auto-enables rerere once `.git/rr-cache` already exists, which
    is exactly the chicken-and-egg case a fresh test repo is in, so a caller
    that forgot to set it first would otherwise get a silently unrecorded
    resolution. The `git config` call is idempotent, so a caller that already
    enabled it (as some do, for their own repo-wide setup) is unaffected.
    """
    run(repo, "config", "rerere.enabled", "true")
    subprocess.run(
        ["git", "merge", "--no-ff", "--no-edit", branch],
        cwd=repo,
        capture_output=True,
        text=True,
        check=False,
    )
    (repo / path).write_text(resolution)
    run(repo, "add", path)
    run(repo, "-c", "commit.gpgsign=false", "commit", "--no-edit")


def branch_touching(repo: Path, branch: str, name: str, content: str) -> None:
    """Create `branch` off master adding a distinct file, then mirror it to `fork/`.

    Mirroring into refs/remotes/fork/* models the real topology: assembly resolves
    features through a remote namespace, never through local branch names.
    """
    run(repo, "checkout", "-q", "-b", branch, "master")
    commit_file(repo, name, content, f"add {name}")
    run(repo, "checkout", "-q", "master")
    run(repo, "update-ref", f"refs/remotes/fork/{branch}", branch)


def conflicting_branch(repo: Path, branch: str) -> None:
    """A branch that edits the same line upstream later edits."""
    run(repo, "checkout", "-q", "-b", branch, "master")
    commit_file(repo, "src.c", "line1\nMINE\n", f"{branch} edits line2")
    run(repo, "checkout", "-q", "master")
    run(repo, "update-ref", f"refs/remotes/fork/{branch}", branch)
    commit_file(repo, "src.c", "line1\nUPSTREAM\n", "upstream edits line2")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """An initialised repo on `master` with a two-line base file."""
    path = tmp_path / "repo"
    path.mkdir()
    run(path, "init", "-q", "-b", "master")
    run(path, "config", "user.email", "test@example.com")
    run(path, "config", "user.name", "Test")
    commit_file(path, "src.c", "line1\nline2\n", "base")
    return path
