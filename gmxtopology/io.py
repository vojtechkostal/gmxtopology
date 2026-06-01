"""Topology file I/O.

`gmxtopology.io` is the place to look for reading and writing topology files.
"""

from itertools import groupby
from operator import attrgetter
from pathlib import Path
from typing import Any, Iterator, Optional

from .topology import Define, MoleculeType, PreprocessorState, RawSection, Topology
from .parser import MOLECULE_SECTIONS, TOPOLOGY_SECTIONS, apply_section_line


SECTION_ALIASES = {
    "dummies3": "virtual_sites3",
}

KNOWN_SECTIONS = {
    "defaults",
    "atomtypes",
    "moleculetype",
    "atoms",
    "exclusions",
    "system",
    "molecules",
    *TOPOLOGY_SECTIONS,
    *MOLECULE_SECTIONS,
}


def _handle_directive(
    line: str,
    fn: Path,
    top: Topology,
    state: PreprocessorState,
    current_section: Optional[str],
    active_mol: Optional[MoleculeType],
) -> tuple[PreprocessorState, Optional[str], Optional[MoleculeType]]:
    directive, *rest = line.split()

    if directive == "#include":
        if len(rest) != 1:
            raise ValueError(
                f"Invalid #include directive in {top.source}: '{line}'"
            )
        inc = rest[0].strip('"<>')
        include = (fn.parent / inc).resolve()
        if state and not include.exists():
            return state, current_section, active_mol
        return read_topology(
            include,
            top,
            initial_state=state,
            initial_section=current_section,
            initial_molecule=active_mol,
            validate_defaults=False,
        )

    if directive == "#define":
        if not rest:
            raise ValueError(
                f"Invalid #define directive in {top.source}: '{line}'. "
                "#define must be followed by a name."
            )
        argument = " ".join(rest[1:]) or None
        define = Define(
            directive=rest[0],
            argument=argument,
            ifdef_state=state,
        )
        top.defines.append(define)
        return state, current_section, active_mol

    if directive in {"#ifdef", "#ifndef"}:
        if len(rest) != 1:
            raise ValueError(
                f"Invalid {directive} directive in {top.source}: '{line}'"
            )
        return (
            (*state, f"{directive.removeprefix('#')} {rest[0]}"),
            current_section,
            active_mol,
        )

    if directive == "#else":
        if rest:
            raise ValueError(f"Invalid #else directive in {top.source}: '{line}'")
        if not state:
            raise ValueError(
                f"#else without matching conditional in {top.source}."
            )
        kind, label = state[-1].split(" ", 1)
        if kind == "else":
            raise ValueError(f"Duplicate #else for '{label}' in {top.source}.")
        return (*state[:-1], f"else {kind} {label}"), current_section, active_mol

    if directive == "#endif":
        if rest:
            raise ValueError(f"Invalid #endif directive in {top.source}: '{line}'")
        if not state:
            raise ValueError(f"#endif without matching conditional in {top.source}.")
        return state[:-1], current_section, active_mol

    raise NotImplementedError(
        f"Unsupported preprocessor directive in {top.source}: '{line}'"
    )


def _directives_for_state(state: str) -> list[str]:
    if state.startswith(("ifdef ", "ifndef ")):
        return [f"#{state}"]
    if state.startswith("else "):
        _, kind, label = state.split(" ", 2)
        return [f"#{kind} {label}", "#else"]
    raise ValueError(f"Unsupported preprocessor state: {state!r}")


def _switch_state(
    lines: list[str],
    active_state: PreprocessorState,
    next_state: PreprocessorState,
) -> PreprocessorState:
    if next_state == active_state:
        return active_state

    if (
        len(active_state) == len(next_state)
        and active_state[:-1] == next_state[:-1]
        and next_state
        and next_state[-1].startswith("else ")
        and next_state[-1] == f"else {active_state[-1]}"
    ):
        lines.append("#else")
        return next_state

    common = 0
    for active, next_ in zip(active_state, next_state):
        if active != next_:
            break
        common += 1

    lines.extend("#endif" for _ in active_state[common:])
    for state in next_state[common:]:
        lines.extend(_directives_for_state(state))
    return next_state


def _iter_section_blocks(
    items: list[Any],
) -> Iterator[tuple[PreprocessorState, str, list[Any]]]:
    if not items:
        return

    groups: dict[int | None, list[Any]] = {}
    for item in items:
        func = item.func if hasattr(item, "func") else None
        groups.setdefault(func, []).append(item)

    for func_items in groups.values():
        for state, state_items in groupby(func_items, key=attrgetter("ifdef_state")):
            block = list(state_items)
            yield state, block[0].header, block


def _write_section_blocks(
    lines: list[str],
    items: list[Any],
    active_state: PreprocessorState,
) -> PreprocessorState:
    for state, header, block_items in _iter_section_blocks(items):
        active_state = _switch_state(lines, active_state, state)
        lines.extend([header, block_items[0].legend])
        lines.extend(str(item) for item in block_items)
        lines.append("")
    return active_state


def _write_raw_section_blocks(
    lines: list[str],
    sections: list[RawSection],
    active_state: PreprocessorState,
) -> PreprocessorState:
    for section in sections:
        active_state = _switch_state(lines, active_state, section.ifdef_state)
        lines.append(section.header)
        lines.extend(section.lines)
        lines.append("")
    return active_state


def _close_state(lines: list[str], active_state: PreprocessorState) -> None:
    lines.extend("#endif" for _ in active_state)


def read_topology(
    fn: Path,
    top: Topology,
    *,
    initial_state: PreprocessorState = (),
    initial_section: Optional[str] = None,
    initial_molecule: Optional[MoleculeType] = None,
    validate_defaults: bool = True,
) -> tuple[PreprocessorState, Optional[str], Optional[MoleculeType]]:
    """
    Reads a topology file and populates the provided Topology object.

    Parameters
    ----------
    fn : Path
        Path to the topology file to be read.
    top : Topology
        An instance of the Topology class to populate with data from the file.

    Notes
    -----
    - The function processes the topology file line by line.
    - It handles conditional blocks (e.g., `#ifdef`, `#ifndef`, `#else`,
      `#endif`) and updates the `ifdef_state` accordingly.
    - Lines starting with `;` are treated as comments and ignored.
    - The function assumes the file is well-formed
      and does not handle all potential edge cases.

    Raises
    ------
    TypeError
        If the input parameters are of incorrect types.
    FileNotFoundError
        If the file specified by `fn` does not exist.
    ValueError
        If the file contains invalid or unexpected content.
    """

    # Validate input types
    if not isinstance(fn, Path):
        raise TypeError("The 'fn' parameter must be of type 'Path'.")
    if not fn.exists():
        raise FileNotFoundError(f"The file '{fn}' does not exist.")
    if not isinstance(top, Topology):
        raise TypeError(
            "The 'top' parameter must be an instance of the 'Topology' class."
        )

    previous_source = top.source
    top.source = fn.resolve()
    current_section = initial_section
    active_mol = initial_molecule
    state = initial_state

    try:
        with fn.open() as f:
            for raw_line in f:
                raw_line = raw_line.rstrip("\n")
                line = raw_line.split(";", 1)[0].strip()

                if not line:
                    continue

                # Preprocessor directives
                if line.startswith("#"):
                    state, current_section, active_mol = _handle_directive(
                        line,
                        fn,
                        top,
                        state,
                        current_section,
                        active_mol,
                    )
                    continue

                # Section header
                if line.startswith("["):
                    section = line.strip("[]").strip()
                    current_section = SECTION_ALIASES.get(section, section)
                    continue

                if not current_section:
                    continue

                if current_section not in KNOWN_SECTIONS:
                    sections = (
                        active_mol.raw_sections
                        if active_mol is not None
                        else top.raw_sections
                    )
                    if (
                        not sections
                        or sections[-1].name != current_section
                        or sections[-1].ifdef_state != state
                    ):
                        sections.append(
                            RawSection(name=current_section, ifdef_state=state)
                        )
                    sections[-1].lines.append(line)
                    continue

                active_mol = apply_section_line(
                    current_section,
                    line.split(),
                    top,
                    active_mol,
                    state,
                )

        if state != initial_state:
            raise ValueError(f"Unclosed conditional block in {top.source}.")
        if validate_defaults and top.defaults is None:
            raise ValueError("No [ defaults ] section found in topology.")
        return state, current_section, active_mol
    finally:
        top.source = previous_source


def write_topology(top: Topology, fn_out: Path, overwrite: bool = False) -> None:
    """Writes the Topology object to a file."""

    fn_out = Path(fn_out)
    if fn_out.exists() and not overwrite:
        raise FileExistsError(
            f"File '{fn_out}' already exists. "
            f"Use overwrite=True to overwrite."
        )

    # Collect lines to write
    lines: list[str] = []

    # --- defines ---
    if top.defines:
        active_state: PreprocessorState = ()
        for define in top.defines:
            active_state = _switch_state(lines, active_state, define.ifdef_state)
            lines.append(str(define))
        _close_state(lines, active_state)
        lines.append("")

    # --- defaults ---
    lines += [top.defaults.header, top.defaults.legend, str(top.defaults), ""]

    used_atomtypes = {
        atomtype.name
        for mol, _ in top.molecules.values()
        for atom in mol.atoms
        for atomtype in (atom.type, atom.type_b)
        if atomtype is not None
    }
    raw_tokens = {
        token.rstrip("\\")
        for section in top.raw_sections
        for line in section.lines
        for token in line.split()
    }
    used_atomtypes.update(
        atomtype.name
        for atomtype in top.atomtypes
        if atomtype.name in raw_tokens
    )

    # --- atomtypes ---
    active_state: PreprocessorState = ()
    atomtypes = [
        atomtype
        for atomtype in top.atomtypes
        if atomtype.name in used_atomtypes
    ]
    active_state = _write_section_blocks(lines, atomtypes, active_state)
    _close_state(lines, active_state)

    # Molecular interactions are written with their resolved parameters.
    # Non-bonded overrides remain global because they cannot be expanded there.
    nonbond_params = [
        param
        for param in top.nonbond_params
        if param.ai.name in used_atomtypes and param.aj.name in used_atomtypes
    ]
    active_state = ()
    active_state = _write_section_blocks(lines, nonbond_params, active_state)
    _close_state(lines, active_state)

    active_state = ()
    active_state = _write_raw_section_blocks(lines, top.raw_sections, active_state)
    _close_state(lines, active_state)

    # --- moleculetypes ---
    for mol, _ in top.molecules.values():

        # moleculetype
        active_state = ()
        active_state = _switch_state(lines, active_state, mol.ifdef_state)
        lines += [mol.header, mol.legend, str(mol), ""]

        active_state = _write_section_blocks(lines, mol.atoms, active_state)

        # interactions
        for section in MOLECULE_SECTIONS:
            params = getattr(mol, section, None)
            if not params:
                continue

            active_state = _write_section_blocks(lines, params, active_state)

        active_state = _write_raw_section_blocks(lines, mol.raw_sections, active_state)

        # exclusions are treated separately
        if mol.exclusions:
            active_state = _write_section_blocks(lines, mol.exclusions, active_state)

        _close_state(lines, active_state)

    # --- system ---
    lines += [top.system.header, top.system.legend, str(top.system), ""]

    # --- molecules ---
    lines += ["[ molecules ]", f"; {'name':>1} {'count':>14}"]
    for mol_name, (_, count) in top.molecules.items():
        lines.append(f"{mol_name:<10} {count:>10}")

    fn_out.write_text("\n".join(lines) + "\n")
