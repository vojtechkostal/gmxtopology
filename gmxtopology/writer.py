"""Write a topology as a single self-contained, flattened ``.top`` file.

Only what the simulated system needs is written: the molecule types listed in
``[ molecules ]``, the atom types they use and the CMAP grids they refer to.
Parameters of molecular interactions were resolved when the topology was read,
so the bonded parameter tables are not written. Conditional blocks and includes
were evaluated and do not appear.
"""

import numbers
from itertools import groupby
from pathlib import Path
from typing import Any, Iterable, Sequence

from .model import Interaction, MoleculeType, Topology, TypeEntry
from .reduce import (
    used_atomtypes,
    used_cmaptypes,
    used_moleculetypes,
    used_nonbond_params,
)

ATOM_LEGEND = ("ai", "aj", "ak", "al", "am")
# Molecular sections in the order they are written.
SECTION_ORDER = (
    "bonds",
    "pairs",
    "pairs_nb",
    "angles",
    "dihedrals",
    "cmap",
    "constraints",
    "settles",
    "exclusions",
    "virtual_sites1",
    "virtual_sites2",
    "virtual_sites3",
    "virtual_sites4",
    "virtual_sitesn",
    "position_restraints",
    "distance_restraints",
    "dihedral_restraints",
    "orientation_restraints",
    "angle_restraints",
    "angle_restraints_z",
)


class _Writer:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.macros: set[str] = set()

    def value(self, value: Any) -> str:
        if isinstance(value, numbers.Integral):
            return str(int(value))
        if isinstance(value, numbers.Real):
            return repr(float(value))
        self.macros.add(str(value).lstrip("-"))
        return str(value)

    def table(
        self, header: str, legend: Sequence[str], rows: Iterable[Sequence[Any]]
    ) -> None:
        """A ``[ header ]`` block with right-aligned columns."""
        cells = [[self.value(v) for v in row] for row in rows]
        if not cells:
            return
        legend = list(legend)[: max(len(row) for row in cells)]
        n_columns = max(len(row) for row in cells)
        widths = [
            max(len(row[i]) for row in [list(legend), *cells] if i < len(row))
            for i in range(n_columns)
        ]

        def align(row: Sequence[str]) -> str:
            return "  ".join(c.rjust(w) for c, w in zip(row, widths)).rstrip()

        self.lines += [
            f"[ {header} ]",
            "; " + align(legend),
            *("  " + align(r) for r in cells),
            "",
        ]

    def interactions(self, section: str, items: list[Interaction]) -> None:
        for _, run in groupby(items, key=lambda item: item.func):
            group = list(run)
            longest = max(group, key=lambda item: len(item.params))
            legend = [*ATOM_LEGEND[: len(longest.atoms)], "func", *longest.params]
            self.table(
                section, legend, ((*i.atoms, i.func, *i.params.values()) for i in group)
            )

    def type_entries(self, section: str, entries: list[TypeEntry]) -> None:
        for _, run in groupby(entries, key=lambda entry: entry.func):
            group = list(run)
            longest = max(group, key=lambda entry: len(entry.params))
            legend = [
                *("type%d" % (i + 1) for i in range(len(longest.types))),
                "func",
                *longest.params,
            ]
            self.table(
                section, legend, ((*e.types, e.func, *e.params.values()) for e in group)
            )

    def molecule(self, mol: MoleculeType) -> None:
        self.table("moleculetype", ["name", "nrexcl"], [(mol.name, mol.nrexcl)])
        rows = []
        for a in mol.atoms:
            row = [a.nr, a.type, a.resnr, a.resname, a.name, a.cgnr, a.charge, a.mass]
            if a.type_b is not None:
                row += [a.type_b, a.charge_b, a.mass_b]
            rows.append(row)
        legend = [
            "nr",
            "type",
            "resnr",
            "residue",
            "atom",
            "cgnr",
            "charge",
            "mass",
            "typeB",
            "chargeB",
            "massB",
        ]
        self.table("atoms", legend, rows)
        for section in SECTION_ORDER:
            if section == "exclusions":
                self.table("exclusions", ["ai", "aj", "..."], mol.exclusions)
            elif mol.interactions.get(section):
                self.interactions(section, mol.interactions[section])


def write_topology(top: Topology, fn_out: Path, overwrite: bool = False) -> None:
    """Write ``top`` to ``fn_out``."""
    fn_out = Path(fn_out)
    if fn_out.exists() and not overwrite:
        raise FileExistsError(
            f"File '{fn_out}' already exists. Use overwrite=True to overwrite."
        )
    if top.defaults is None:
        raise ValueError("Cannot write a topology without [ defaults ].")

    w = _Writer()
    d = top.defaults
    legend = ["nbfunc", "comb-rule", "gen-pairs", "fudgeLJ", "fudgeQQ"]
    row = [d.nbfunc, d.comb_rule, d.gen_pairs, d.fudgeLJ, d.fudgeQQ]
    if d.n is not None:
        legend, row = [*legend, "N"], [*row, d.n]
    w.table("defaults", legend, [row])

    atomtypes = used_atomtypes(top)
    rows, legends = [], []
    for t in atomtypes:
        optional = [("bonded_type", t.bonded_type), ("at.num", t.atnum)]
        names = [name for name, value in optional if value is not None]
        values = [value for _, value in optional if value is not None]
        legends.append(["name", *names, "mass", "charge", "ptype", "sigma", "epsilon"])
        rows.append([t.name, *values, t.mass, t.charge, t.ptype, t.sigma, t.epsilon])
    w.table("atomtypes", max(legends, key=len, default=[]), rows)

    nonbond = used_nonbond_params(top, atomtypes)
    w.type_entries("nonbond_params", nonbond)

    cmaptypes = used_cmaptypes(top)
    if cmaptypes:
        w.lines.append("[ cmaptypes ]")
        for c in cmaptypes:
            values = [repr(float(v)) for v in c.values]
            w.lines.append(
                " ".join([*c.types, str(c.func), str(c.nx), str(c.ny)]) + " \\"
            )
            chunks = [values[i : i + 10] for i in range(0, len(values), 10)]
            w.lines += [
                " ".join(chunk) + (" \\" if i < len(chunks) - 1 else "")
                for i, chunk in enumerate(chunks)
            ]
        w.lines.append("")

    for mol in used_moleculetypes(top):
        w.molecule(mol)

    w.lines += ["[ system ]", top.system or "system", ""]
    w.table(
        "molecules", ["name", "count"], [(m.name, count) for m, count in top.molecules]
    )

    if any(top.intermolecular.values()):
        w.lines += ["[ intermolecular_interactions ]"]
        for section, items in top.intermolecular.items():
            w.interactions(section, items)

    # Macros that parameters refer to must be defined in the output.
    defines = {d.directive: d for d in top.defines}
    defines.update(
        {n: d for n, d in top.external_defines.items() if d.argument is not None}
    )
    header = [str(defines[name]) for name in sorted(w.macros) if name in defines]
    text = "\n".join([*header, *([""] if header else []), *w.lines])
    fn_out.write_text(text.rstrip("\n") + "\n")
