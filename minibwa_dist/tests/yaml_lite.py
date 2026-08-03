"""A minimal, stdlib-only YAML subset parser for the GitHub Actions workflows.

Not a general YAML implementation -- deliberately so. It exists for exactly one
job: turn a workflow file into nested `dict`/`list`/`str` so `test_workflows.py`
can assert against `jobs.<job>.steps[N].run` (or `.name`, `.if`, ...) instead of
scanning the raw file text. That distinction is the whole point: a real parser
drops `#` comments between YAML nodes for free, and confines a match to the one
step it names instead of "matched somewhere in the file". `pyyaml` is not a
dependency of this project and is not available through the pinned pixi
environment these tests run under, so this covers the subset the workflows
actually use: block mappings and sequences, plain/quoted scalars, and literal
(`|`) or folded (`>`) block scalars with optional chomping indicators. No flow
collections, anchors, aliases, or multi-document streams -- none of the
workflows use them, and this parser does not try to guess if one shows up.

All mapping keys and plain scalar values are kept as plain `str`. Real YAML 1.1
would coerce a bare `on:` key to the boolean `True` (the classic GitHub Actions
gotcha), and bare `true`/`false` values to `bool` -- neither transformation is
wanted here, so it is simply never done.
"""

from __future__ import annotations

import re

_BLOCK_SCALAR_STYLES = frozenset({"|", "|-", "|+", ">", ">-", ">+"})
_TRAILING_COMMENT = re.compile(r"(?<=\s)#.*$")


class YamlLiteError(ValueError):
    """The input used a construct this parser does not understand."""


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _strip_trailing_comment(value: str) -> str:
    """Drop a ` # ...` suffix from a plain scalar, e.g. `uses: foo@sha # v7.0.1`.

    Only ever applied to plain (unquoted, non-block) scalars: none of those
    values in these workflows legitimately contain a `#`.
    """
    return _TRAILING_COMMENT.sub("", value).rstrip()


def _unquote(value: str) -> str:
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1].replace('\\"', '"')
    if len(value) >= 2 and value[0] == "'" and value[-1] == "'":
        return value[1:-1].replace("''", "'")
    return _strip_trailing_comment(value)


class _Parser:
    """Single-pass recursive-descent parser over the document's raw lines."""

    def __init__(self, text: str) -> None:
        self._lines = text.splitlines()
        self._idx = 0

    def parse(self) -> dict[str, object] | list[object] | None:
        return self._parse_block(min_indent=0)

    # -- line cursor -----------------------------------------------------

    def _peek_real_line(self) -> str | None:
        """Advance past blank lines and whole-line YAML comments; return the
        next substantive line without consuming it. `None` at end of input.

        Never called while collecting a block scalar -- that path reads
        `self._lines` directly, because blank lines and `#`-led lines inside a
        `run: |` block are shell content, not YAML comments.
        """
        while self._idx < len(self._lines):
            line = self._lines[self._idx]
            stripped = line.strip()
            if stripped == "" or stripped.startswith("#"):
                self._idx += 1
                continue
            return line
        return None

    # -- block dispatch ----------------------------------------------------

    def _parse_block(self, min_indent: int) -> dict[str, object] | list[object] | None:
        line = self._peek_real_line()
        if line is None:
            return None
        indent = _indent(line)
        if indent < min_indent:
            return None
        content = line[indent:]
        if content == "-" or content.startswith("- "):
            return self._parse_sequence(indent)
        return self._parse_mapping(indent)

    def _parse_sequence(self, indent: int) -> list[object]:
        result: list[object] = []
        while True:
            line = self._peek_real_line()
            if line is None or _indent(line) != indent:
                break
            content = line[indent:]
            if not (content == "-" or content.startswith("- ")):
                break
            after_dash = content[1:]
            if after_dash.strip() == "":
                self._idx += 1
                result.append(self._parse_block(min_indent=indent + 1))
                continue
            leading = len(after_dash) - len(after_dash.lstrip(" "))
            child_indent = indent + 1 + leading
            # Rewrite "- key: value" in place into a plain "key: value" line at
            # `child_indent`, so the mapping parser below can consume it (and
            # any siblings already indented to match it) uniformly. The cursor
            # is deliberately not advanced -- the rewritten line still needs
            # to be read as the first line of that mapping.
            self._lines[self._idx] = " " * child_indent + after_dash.lstrip(" ")
            result.append(self._parse_mapping(child_indent))
        return result

    def _parse_mapping(self, indent: int) -> dict[str, object]:
        result: dict[str, object] = {}
        while True:
            line = self._peek_real_line()
            if line is None or _indent(line) != indent:
                break
            content = line[indent:]
            if content.startswith("- "):
                break
            key, sep, rest = content.partition(":")
            if not sep:
                raise YamlLiteError(f"expected 'key: value', got: {line!r}")
            key = key.strip()
            self._idx += 1
            value_str = rest.strip()
            if value_str == "":
                result[key] = self._parse_block(min_indent=indent + 1)
            elif value_str in _BLOCK_SCALAR_STYLES:
                result[key] = self._parse_block_scalar(indent, chomp=value_str[1:])
            else:
                result[key] = _unquote(value_str)
        return result

    def _parse_block_scalar(self, parent_indent: int, chomp: str) -> str:
        """Collect a literal/folded block scalar's raw lines, dedented.

        Folding (`>`) is not actually performed -- nothing this suite asserts
        on reads a folded value's content, only a literal `run:` script's, so
        treating both styles as "keep the lines verbatim" is a deliberate
        simplification, not an oversight. Chomping is likewise approximated:
        `-`/`+` distinguish only by how many trailing blank lines survive, and
        no assertion here depends on that either, so trailing blanks are
        always trimmed.
        """
        block_indent: int | None = None
        collected: list[str] = []
        while self._idx < len(self._lines):
            raw = self._lines[self._idx]
            if raw.strip() == "":
                collected.append("")
                self._idx += 1
                continue
            cur_indent = _indent(raw)
            if block_indent is None:
                if cur_indent <= parent_indent:
                    break
                block_indent = cur_indent
            if cur_indent < block_indent:
                break
            collected.append(raw[block_indent:])
            self._idx += 1
        while collected and collected[-1] == "":
            collected.pop()
        return "\n".join(collected)


def parse(text: str) -> dict[str, object]:
    """Parse a GitHub Actions workflow file's text into nested dict/list/str."""
    result = _Parser(text).parse()
    if not isinstance(result, dict):
        raise YamlLiteError(f"expected a top-level mapping, got {type(result).__name__}")
    return result


def step(workflow: dict[str, object], job: str, *, name: str) -> dict[str, object]:
    """The `jobs.<job>.steps[]` entry whose `name` field equals `name`.

    Raises with the available names on a miss, rather than returning `None` --
    a step lookup that silently finds nothing is exactly the "matched
    somewhere else, or nowhere" failure mode this module exists to rule out.
    """
    jobs = workflow.get("jobs")
    if not isinstance(jobs, dict) or job not in jobs:
        raise YamlLiteError(f"no job '{job}'; jobs present: {sorted(jobs) if jobs else []}")
    steps = jobs[job].get("steps")
    if not isinstance(steps, list):
        raise YamlLiteError(f"job '{job}' has no steps list")
    for entry in steps:
        if isinstance(entry, dict) and entry.get("name") == name:
            return entry
    names = [entry.get("name") for entry in steps if isinstance(entry, dict)]
    raise YamlLiteError(f"no step named {name!r} in job '{job}'; steps: {names}")
