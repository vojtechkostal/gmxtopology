"""A small C-like preprocessor for GROMACS topology files.

It supports the directives of the GROMACS preprocessor (``#include``,
``#define``, ``#undef``, ``#ifdef``, ``#ifndef``, ``#else``, ``#endif``),
strips ``;`` comments and joins ``\\`` continuation lines. Conditionals are
evaluated against the macros that are defined at that point, so only the
active branch is yielded.
"""

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Sequence

MAX_INCLUDE_DEPTH = 64


class TopologyError(ValueError):
    """An error in a topology file, with its location when known."""


@dataclass(frozen=True)
class Line:
    """One logical topology line and where it came from."""

    text: str
    file: Path
    lineno: int

    def __str__(self) -> str:
        return f"{self.file}:{self.lineno}"


@dataclass
class Define:
    """A ``#define`` macro. ``argument`` is the replacement text, if any."""

    directive: str
    argument: str | int | float | None = None

    def update(self, *, argument: str | int | float | None) -> None:
        self.argument = argument

    def __str__(self) -> str:
        if self.argument is None:
            return f"#define {self.directive}"
        return f"#define {self.directive} {self.argument}"


def gromacs_include_dirs() -> list[Path]:
    """Directories GROMACS searches for force-field includes."""
    dirs = [Path(p) for p in os.environ.get("GMXLIB", "").split(os.pathsep) if p]
    if os.environ.get("GMXDATA"):
        dirs.append(Path(os.environ["GMXDATA"]) / "top")
    gmx = shutil.which("gmx")
    if gmx:
        dirs.append(Path(gmx).resolve().parents[1] / "share" / "gromacs" / "top")
    return dirs


def _find_include(name: str, current: Path, search: Sequence[Path]) -> Path | None:
    for directory in (current.parent, *search):
        candidate = directory / name
        if candidate.is_file():
            return candidate.resolve()
    return None


def _logical_lines(path: Path) -> Iterator[tuple[int, str]]:
    """Comment-stripped lines with ``\\`` continuations joined."""
    pending, start = "", 0
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        text = raw.split(";", 1)[0].rstrip()
        if not pending:
            start = lineno
        if text.endswith("\\"):
            pending += text[:-1] + " "
            continue
        text, pending = pending + text, ""
        if text.strip():
            yield start, text.strip()
    if pending.strip():
        yield start, pending.strip()


def preprocess(
    path: Path,
    macros: dict[str, Define],
    defined: list[Define],
    include_dirs: Sequence[Path] = (),
    _depth: int = 0,
) -> Iterator[Line]:
    """Yield the active logical lines of ``path`` and everything it includes.

    ``macros`` holds the macros currently in force and is updated in place.
    Macros defined by the files themselves are also appended to ``defined``
    (replacing an earlier definition of the same name).
    """
    if _depth > MAX_INCLUDE_DEPTH:
        raise TopologyError(f"Include nesting is too deep at {path}.")
    if not path.is_file():
        raise FileNotFoundError(f"Topology file '{path}' does not exist.")

    # One entry per open conditional: (this branch is active, a branch was taken).
    stack: list[tuple[bool, bool]] = []

    for lineno, text in _logical_lines(path):
        line = Line(text, path, lineno)
        active = all(branch for branch, _ in stack)

        if not text.startswith("#"):
            if active:
                yield line
            continue

        directive, *rest = text.split()
        if directive in ("#ifdef", "#ifndef"):
            if len(rest) != 1:
                raise TopologyError(f"Invalid {directive} at {line}: '{text}'.")
            taken = (rest[0] in macros) == (directive == "#ifdef")
            stack.append((taken, taken))
        elif directive == "#else":
            if not stack or rest:
                raise TopologyError(f"Invalid or unmatched #else at {line}.")
            _, taken = stack[-1]
            stack[-1] = (not taken, True)
        elif directive == "#endif":
            if not stack or rest:
                raise TopologyError(f"Invalid or unmatched #endif at {line}.")
            stack.pop()
        elif not active:
            continue
        elif directive == "#define":
            if not rest:
                raise TopologyError(f"#define without a name at {line}.")
            define = Define(rest[0], " ".join(rest[1:]) or None)
            macros[define.directive] = define
            defined[:] = [d for d in defined if d.directive != define.directive]
            defined.append(define)
        elif directive == "#undef":
            if len(rest) != 1:
                raise TopologyError(f"Invalid #undef at {line}: '{text}'.")
            macros.pop(rest[0], None)
            defined[:] = [d for d in defined if d.directive != rest[0]]
        elif directive == "#include":
            if len(rest) != 1:
                raise TopologyError(f"Invalid #include at {line}: '{text}'.")
            name = rest[0].strip('"<>')
            include = _find_include(
                name, path, [*include_dirs, *gromacs_include_dirs()]
            )
            if include is None:
                raise FileNotFoundError(f"Cannot find include '{name}' at {line}.")
            yield from preprocess(include, macros, defined, include_dirs, _depth + 1)
        else:
            raise TopologyError(
                f"Unsupported preprocessor directive at {line}: '{text}'."
            )

    if stack:
        raise TopologyError(f"Unclosed conditional block in {path}.")


def initial_macros(defines: Mapping[str, object] | None) -> dict[str, Define]:
    """Macros supplied by the caller, like ``-DNAME`` or ``-DNAME=value``."""
    return {
        name: Define(name, None if value in (None, True, "") else value)
        for name, value in (defines or {}).items()
    }
