"""Topology file I/O.

`gmxtopology.io` is the place to look for reading and writing topology files.
"""

from pathlib import Path
from typing import Optional, Dict, List, Any

from .topology import Topology, Define, MoleculeType
from .parser import MOLECULE_SECTIONS, apply_section_line


SECTION_ALIASES = {
    "dummies3": "virtual_sites3",
}


def _normalize_section(section: Optional[str]) -> Optional[str]:
    if section is None:
        return None
    return SECTION_ALIASES.get(section, section)


def _handle_directive(
    line: str,
    fn: Path,
    top: Topology,
    ifdef_block: Optional[str],
    ifdef_state: str,
) -> tuple[Optional[str], str]:
    directive, *rest = line.split()

    if directive == "#include":
        inc = rest[0].strip('"')
        read_topology((fn.parent / inc).resolve(), top)
        return ifdef_block, ifdef_state

    if directive == "#define":
        if len(rest) != 2:
            raise ValueError(
                f"Invalid #define directive in {top.source}: '{line}'"
                "#define must be followed by name and value."
            )
        define = Define(directive=rest[0], argument=rest[1])
        top.defines.add(define)
        return ifdef_block, ifdef_state

    if directive == "#ifdef":
        if len(rest) != 1:
            raise ValueError(
                f"Invalid #ifdef directive in {top.source}: '{line}'"
            )
        ifdef_block = rest[0]
        return ifdef_block, f"ifdef {ifdef_block}"

    if directive == "#else":
        if ifdef_block is None:
            raise ValueError(
                f"#else without matching #ifdef in {top.source}."
            )
        return ifdef_block, f"else {ifdef_block}"

    if directive == "#endif":
        return None, "free"

    return ifdef_block, ifdef_state


def _group_by_func(items: List) -> Dict[int, List]:
    """Group items by their 'func' attribute."""
    groups: Dict[int, List] = {}
    for item in items:
        groups.setdefault(item.func, []).append(item)
    return groups


def _directive_for_state(state: str) -> str:
    if state.startswith("ifdef "):
        return f"#{state}"
    if state.startswith("else "):
        return "#else"
    raise ValueError(f"Unsupported preprocessor state: {state!r}")


def _switch_state(
    lines: list[str],
    active_state: str | None,
    next_state: str,
) -> str | None:
    if next_state == active_state:
        return active_state

    if active_state is None:
        if next_state != "free":
            lines.append(_directive_for_state(next_state))
            return next_state
        return None

    if next_state == "free":
        lines.append("#endif")
        return None

    active_kind, active_label = active_state.split(" ", 1)
    next_kind, next_label = next_state.split(" ", 1)
    if active_kind == "ifdef" and next_kind == "else" and active_label == next_label:
        lines.append("#else")
        return next_state

    lines.append("#endif")
    lines.append(_directive_for_state(next_state))
    return next_state


def _partition_by_state(items: list[Any]) -> list[list[Any]]:
    groups: list[list[Any]] = []
    for item in items:
        if not groups or groups[-1][0].ifdef_state != item.ifdef_state:
            groups.append([item])
        else:
            groups[-1].append(item)
    return groups


def _iter_section_blocks(items: list[Any]) -> list[tuple[str, str, list[Any]]]:
    blocks: list[tuple[str, str, list[Any]]] = []
    if not items:
        return blocks

    if hasattr(items[0], "func"):
        for func_items in _group_by_func(items).values():
            for state_items in _partition_by_state(func_items):
                blocks.append(
                    (state_items[0].ifdef_state, state_items[0].header, state_items)
                )
        return blocks

    for state_items in _partition_by_state(items):
        blocks.append((state_items[0].ifdef_state, state_items[0].header, state_items))
    return blocks


def _write_section_blocks(
    lines: list[str],
    items: list[Any],
    active_state: str | None,
) -> str | None:
    for state, header, block_items in _iter_section_blocks(items):
        active_state = _switch_state(lines, active_state, state)
        lines.extend([header, block_items[0].legend])
        lines.extend(str(item) for item in block_items)
        lines.append("")
    return active_state


def _close_state(lines: list[str], active_state: str | None) -> None:
    if active_state is not None:
        lines.append("#endif")


def read_topology(fn: Path, top: Topology) -> None:
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
    - It handles conditional blocks (e.g., `#ifdef`, `#else`, `#endif`)
      and updates the `ifdef_state` accordingly.
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

    top.source = fn.resolve()
    current_section: Optional[str] = None
    active_mol: Optional[MoleculeType] = None
    ifdef_block: Optional[str] = None
    ifdef_state = "free"

    with fn.open() as f:
        for raw_line in f:
            raw_line = raw_line.rstrip("\n")
            line = raw_line.split(";", 1)[0].strip()

            if not line:
                continue

            # Preprocessor directives
            if line.startswith("#"):
                ifdef_block, ifdef_state = _handle_directive(
                    line, fn, top, ifdef_block, ifdef_state
                )
                continue

            # Section header
            if line.startswith("["):
                current_section = _normalize_section(line.strip("[]").strip())
                continue

            if not current_section:
                continue

            active_mol = apply_section_line(
                current_section,
                line.split(),
                top,
                active_mol,
                ifdef_state,
            )

        if top.defaults is None:
            raise ValueError("No [ defaults ] section found in topology.")


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
        for define in top.defines:
            lines.append(str(define))
        lines.append("")

    # --- defaults ---
    lines += [top.defaults.header, top.defaults.legend, str(top.defaults), ""]

    # --- atomtypes ---
    used_atomtypes = list(dict.fromkeys(
        atom.type
        for mol, _ in top.molecules.values()
        for atom in mol.atoms
    ))
    active_state: str | None = None
    active_state = _write_section_blocks(lines, used_atomtypes, active_state)
    _close_state(lines, active_state)

    # --- nonbonded ---
    nonbond_params = getattr(top, "nonbond_params", [])
    if nonbond_params:
        used_nonbond_params = list({
            nb
            for nb in nonbond_params
            if nb.ai in used_atomtypes and nb.aj in used_atomtypes
        })
        active_state = None
        active_state = _write_section_blocks(lines, used_nonbond_params, active_state)
        _close_state(lines, active_state)

    # --- moleculetypes ---
    for mol, _ in top.molecules.values():

        # moleculetype
        lines += [mol.header, mol.legend, str(mol), ""]

        active_state = None
        active_state = _write_section_blocks(lines, mol.atoms, active_state)

        # interactions
        for section in MOLECULE_SECTIONS:
            params = getattr(mol, section, None)
            if not params:
                continue

            active_state = _write_section_blocks(lines, params, active_state)

        # exclusions are treated separately
        if mol.exclusions:
            active_state = _write_section_blocks(lines, mol.exclusions, active_state)

        _close_state(lines, active_state)

    # --- system ---
    lines += [top.system.header, top.system.legend, str(top.system), ""]

    # --- molecules ---
    lines += ["[ molecules ]", f"; {'name':>1} {'count':>14}"]
    for mol_name, (mol, count) in top.molecules.items():
        lines.append(f"{mol_name:<10} {count:>10}")

    fn_out.write_text("\n".join(lines) + "\n")
