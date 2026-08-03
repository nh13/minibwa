"""A thin, testable wrapper around the git CLI."""

from __future__ import annotations

import subprocess
from pathlib import Path


class GitError(subprocess.CalledProcessError):
    """A git command that failed, carrying git's own explanation.

    `subprocess.CalledProcessError` renders as the argv and the exit status and
    nothing else. With `capture_output=True` the reason git actually gave is
    captured into `stderr` and then never shown, so an unattended run files an
    issue pointing at a log that names a command but not a cause. Subclassing
    keeps every `except subprocess.CalledProcessError` working.
    """

    def __str__(self) -> str:
        detail = (self.stderr or "").strip()
        return f"{super().__str__()}\n{detail}" if detail else super().__str__()


def git(repo: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    """Run `git *args` in `repo`, capturing output as text."""
    result = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise GitError(result.returncode, result.args, result.stdout, result.stderr)
    return result


def git_ok(repo: Path, *args: str) -> bool:
    """True when `git *args` exits zero."""
    return git(repo, *args, check=False).returncode == 0


def stdout(repo: Path, *args: str) -> str:
    """Stripped stdout of a successful git command."""
    return git(repo, *args).stdout.strip()
