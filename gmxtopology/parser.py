"""Parse preprocessed topology lines into a :class:`Topology`."""

import warnings
from pathlib import Path
from typing import Any, Iterable

from .lookup import TypeIndex
from .model import (
    Atom,
    AtomType,
    CmapType,
    Defaults,
    Interaction,
    MoleculeType,
    Topology,
    TypeEntry,
)
from .preprocess import Line, TopologyError, initial_macros, preprocess
from .specs import (
    ALIASES,
    DEFAULT_FUNC_SECTIONS,
    INTERMOLECULAR_SECTIONS,
    LOOKUP_SECTIONS,
    MOLECULE_SECTIONS,
    SPECS,
    TYPE_SECTIONS,
    parse_params,
    parse_value,
)

# Parameter tables and the molecular directive they supply defaults for.
TABLE_FOR = {section: table for table, section in TYPE_SECTIONS.items()}
IGNORED_SECTIONS = {"implicit_genborn_params", "implicit_surface_params"}
PTYPES = {"A", "S", "V", "D"}


def read_topology(
    top: Topology, defines: dict[str, Any] | None, include_dirs: list[Path]
) -> None:
    """Read ``top.source`` and everything it includes into ``top``."""
    macros = initial_macros(defines)
    top.external_defines = dict(macros)
    lines = preprocess(top.source, macros, top.defines, include_dirs)
    _Reader(top, macros).read(lines)


class _Reader:
    def __init__(self, top: Topology, macros: dict) -> None:
        self.top = top
        self.macros = macros
        self.section: str | None = None
        self.mol: MoleculeType | None = None
        self.intermolecular = False
        self._index: TypeIndex | None = None
        self._atomtypes_by_name: dict[str, AtomType] = {}

    def read(self, lines: Iterable[Line]) -> None:
        for line in lines:
            if line.text.startswith("["):
                self._enter_section(line)
            elif self.section is not None:  # text before any section is ignored
                self._data(line, self._expand(line.text.split()))

        if self.top.defaults is None:
            raise TopologyError(f"No [ defaults ] section found in {self.top.source}.")

    # -- sections ---------------------------------------------------------

    def _enter_section(self, line: Line) -> None:
        name = line.text.strip("[]").strip()
        name = ALIASES.get(name, name)
        if name in IGNORED_SECTIONS:
            warnings.warn(f"Ignoring obsolete [ {name} ] at {line}.", stacklevel=2)
        if name == "intermolecular_interactions":
            self.intermolecular = True
        self.section = name

    def _data(self, line: Line, tokens: list[str]) -> None:
        section = self.section
        if section in IGNORED_SECTIONS:
            return
        elif section == "intermolecular_interactions":
            raise TopologyError(f"Expected a directive after {section} at {line}.")
        elif self.intermolecular:
            self._interaction(line, tokens)
        elif section in _HANDLERS:
            getattr(self, _HANDLERS[section])(line, tokens)
        elif section in SPECS and section in MOLECULE_SECTIONS:
            self._interaction(line, tokens)
        elif section in TYPE_SECTIONS:
            self._type_entry(line, tokens)
        else:
            raise TopologyError(f"Unknown directive [ {section} ] at {line}.")

    def _expand(self, tokens: list[str]) -> list[str]:
        """Replace macros that stand for several values by those values."""
        expanded: list[str] = []
        for token in tokens:
            macro = self.macros.get(token)
            values = str(macro.argument).split() if macro and macro.argument else []
            expanded.extend(values if len(values) > 1 else [token])
        return expanded

    def _int(self, token: str, line: Line, what: str) -> int:
        try:
            return int(token)
        except ValueError:
            raise TopologyError(f"Invalid {what} '{token}' at {line}.") from None

    def _float(self, token: str, line: Line) -> float | str:
        return parse_value(token, float, self.macros, str(line))

    def _molecule(self, line: Line) -> MoleculeType:
        if self.mol is None:
            raise TopologyError(
                f"[ {self.section} ] before [ moleculetype ] at {line}."
            )
        return self.mol

    def _global_only(self, line: Line) -> None:
        if self.mol is not None:
            raise TopologyError(
                f"[ {self.section} ] after [ moleculetype ] {self.mol.name} at {line}; "
                "parameter tables must come first."
            )

    # -- global sections --------------------------------------------------

    def _defaults(self, line: Line, tokens: list[str]) -> None:
        self._global_only(line)
        if self.top.defaults is not None:
            raise TopologyError(f"Multiple [ defaults ] sections at {line}.")
        if not 2 <= len(tokens) <= 6:
            raise TopologyError(f"[ defaults ] expects 2 to 6 values at {line}.")
        d = Defaults(
            nbfunc=self._int(tokens[0], line, "nbfunc"),
            comb_rule=self._int(tokens[1], line, "comb-rule"),
        )
        if len(tokens) > 2:
            d.gen_pairs = tokens[2]
        if len(tokens) > 3:
            d.fudgeLJ = self._float(tokens[3], line)
        if len(tokens) > 4:
            d.fudgeQQ = self._float(tokens[4], line)
        if len(tokens) > 5:
            d.n = self._int(tokens[5], line, "N")
        if d.nbfunc not in (1, 2) or d.comb_rule not in (1, 2, 3):
            raise TopologyError(f"Invalid nbfunc/comb-rule in [ defaults ] at {line}.")
        if d.gen_pairs not in ("yes", "no"):
            raise TopologyError(f"gen-pairs must be 'yes' or 'no' at {line}.")
        self.top.defaults = d

    def _atomtypes(self, line: Line, tokens: list[str]) -> None:
        self._global_only(line)
        # name [bonded_type] [atnum] mass charge ptype sigma epsilon: locate
        # the one-letter particle type to find which optional columns exist.
        ptype_at = next(
            (i for i in (3, 4, 5) if len(tokens) == i + 3 and tokens[i] in PTYPES),
            None,
        )
        if ptype_at is None:
            raise TopologyError(f"Cannot interpret [ atomtypes ] line at {line}.")
        extra = tokens[1 : ptype_at - 2]
        bonded_type = atnum = None
        if len(extra) == 2:
            bonded_type, atnum = extra[0], self._int(extra[1], line, "atomic number")
        elif extra and extra[0].isdigit():
            atnum = int(extra[0])
        elif extra:
            bonded_type = extra[0]
        mass, charge, _, sigma, epsilon = tokens[ptype_at - 2 :]
        atomtype = AtomType(
            name=tokens[0],
            bonded_type=bonded_type,
            atnum=atnum,
            mass=self._float(mass, line),
            charge=self._float(charge, line),
            ptype=tokens[ptype_at],
            sigma=self._float(sigma, line),
            epsilon=self._float(epsilon, line),
        )
        old = self._atomtypes_by_name.get(atomtype.name)
        if old is None:
            self.top.atomtypes.append(atomtype)
        else:
            warnings.warn(
                f"Atom type {atomtype.name} redefined at {line}.", stacklevel=2
            )
            self.top.atomtypes[self.top.atomtypes.index(old)] = atomtype
        self._atomtypes_by_name[atomtype.name] = atomtype

    def _type_entry(self, line: Line, tokens: list[str]) -> None:
        self._global_only(line)
        section = self.section
        n = SPECS[TYPE_SECTIONS[section]].n_atoms
        if (
            section == "dihedraltypes"
            and len(tokens) > 2
            and tokens[2].isdigit()
            and len(tokens[2]) == 1
        ):
            # Two types: the central pair (the outer pair for impropers).
            func = self._int(tokens[2], line, "function type")
            first, last = tokens[:2]
            types = (
                (first, "X", "X", last) if func in (2, 4) else ("X", first, last, "X")
            )
            rest = tokens[3:]
        else:
            if len(tokens) <= n:
                raise TopologyError(f"Too few values in [ {section} ] at {line}.")
            types, func, rest = (
                tuple(tokens[:n]),
                self._int(tokens[n], line, "function type"),
                tokens[n + 1 :],
            )
        params = parse_params(
            TYPE_SECTIONS[section], func, rest, self.macros, str(line)
        )
        self.top.types.setdefault(section, []).append(
            TypeEntry(section, types, func, params)
        )
        self._index = None

    def _cmaptypes(self, line: Line, tokens: list[str]) -> None:
        self._global_only(line)
        if len(tokens) < 8:
            raise TopologyError(f"Too few values in [ cmaptypes ] at {line}.")
        func, nx, ny = (self._int(t, line, "CMAP value") for t in tokens[5:8])
        values = [self._float(t, line) for t in tokens[8:]]
        if len(values) != nx * ny:
            raise TopologyError(
                f"[ cmaptypes ] at {line} has {len(values)} values, expected {nx * ny}."
            )
        self.top.cmaptypes.append(CmapType(tuple(tokens[:5]), func, nx, ny, values))

    # -- molecules --------------------------------------------------------

    def _moleculetype(self, line: Line, tokens: list[str]) -> None:
        if len(tokens) != 2:
            raise TopologyError(
                f"[ moleculetype ] expects a name and nrexcl at {line}."
            )
        name = tokens[0]
        if any(m.name.casefold() == name.casefold() for m in self.top.moleculetypes):
            raise TopologyError(f"Duplicate moleculetype '{name}' at {line}.")
        self.mol = MoleculeType(name, self._int(tokens[1], line, "nrexcl"))
        self.top.moleculetypes.append(self.mol)

    def _atoms(self, line: Line, tokens: list[str]) -> None:
        mol = self._molecule(line)
        if len(tokens) not in (6, 7, 8, 10, 11):
            raise TopologyError(f"Unexpected number of columns in [ atoms ] at {line}.")
        nr = self._int(tokens[0], line, "atom number")
        if nr != len(mol.atoms) + 1:
            raise TopologyError(
                f"Atom number {nr} in molecule {mol.name} at {line} is not "
                f"{len(mol.atoms) + 1}; atoms must be numbered consecutively."
            )
        atomtype = self._atomtype(tokens[1], line)
        b = tokens[-3:] if len(tokens) in (10, 11) else [None] * 3
        if b[0] is not None:
            self._atomtype(b[0], line)
        mol.atoms.append(
            Atom(
                nr=nr,
                type=tokens[1],
                resnr=self._int(tokens[2], line, "residue number"),
                resname=tokens[3],
                name=tokens[4],
                cgnr=self._int(tokens[5], line, "charge group"),
                charge=self._float(tokens[6], line)
                if len(tokens) > 6
                else atomtype.charge,
                mass=self._float(tokens[7], line)
                if len(tokens) in (8, 11)
                else atomtype.mass,
                type_b=b[0],
                charge_b=self._float(b[1], line) if b[1] else None,
                mass_b=self._float(b[2], line) if b[2] else None,
            )
        )

    def _atomtype(self, name: str, line: Line) -> AtomType:
        try:
            return self._atomtypes_by_name[name]
        except KeyError:
            raise TopologyError(f"Unknown atom type '{name}' at {line}.") from None

    def _exclusions(self, line: Line, tokens: list[str]) -> None:
        mol = self._molecule(line)
        mol.exclusions.append(tuple(self._atom_numbers(tokens, line, mol)))

    def _atom_numbers(self, tokens: list[str], line: Line, mol: MoleculeType | None):
        numbers = [self._int(token, line, "atom number") for token in tokens]
        if mol is not None and not all(1 <= nr <= len(mol.atoms) for nr in numbers):
            raise TopologyError(
                f"Atom number outside 1..{len(mol.atoms)} of molecule {mol.name} at {line}."
            )
        return numbers

    def _interaction(self, line: Line, tokens: list[str]) -> None:
        section = self.section
        spec = SPECS.get(section)
        if spec is None or section == "nonbond_params":
            raise TopologyError(f"Unknown directive [ {section} ] at {line}.")
        if self.intermolecular and section not in INTERMOLECULAR_SECTIONS:
            raise TopologyError(
                f"[ {section} ] is not allowed in intermolecular interactions at {line}."
            )
        mol = None if self.intermolecular else self._molecule(line)
        n = spec.n_atoms
        if len(tokens) < n:
            raise TopologyError(
                f"Expected {n} atom numbers in [ {section} ] at {line}."
            )
        atoms = tuple(self._atom_numbers(tokens[:n], line, mol))
        if len(tokens) == n and section in DEFAULT_FUNC_SECTIONS:
            func, rest = 1, []
        elif len(tokens) == n:
            raise TopologyError(f"Missing function type in [ {section} ] at {line}.")
        else:
            func, rest = self._int(tokens[n], line, "function type"), tokens[n + 1 :]

        if rest:
            params = parse_params(section, func, rest, self.macros, str(line))
            found = [Interaction(section, atoms, func, params)]
        elif (
            func in spec.funcs
            and not spec.funcs[func].names
            and not spec.funcs[func].rest
        ):
            found = [Interaction(section, atoms, func)]
        elif section in LOOKUP_SECTIONS and mol is not None:
            found = self._default_parameters(line, mol, section, atoms, func)
        else:
            raise TopologyError(f"Parameters are required for [ {section} ] at {line}.")

        target = self.top.intermolecular if self.intermolecular else mol.interactions
        target.setdefault(section, []).extend(found)

    def _default_parameters(
        self,
        line: Line,
        mol: MoleculeType,
        section: str,
        atoms: tuple[int, ...],
        func: int,
    ) -> list[Interaction]:
        if self._index is None:
            self._index = TypeIndex(self.top.types)
        table = TABLE_FOR[section]
        records = [mol.atoms[nr - 1] for nr in atoms]

        def names(types: list[str]) -> tuple[str, ...]:
            # Pair parameters are tabulated by atom type, all others by bonded type.
            if section == "pairs":
                return tuple(types)
            return tuple(self._atomtypes_by_name[t].btype for t in types)

        found = self._index.find(table, names([a.type for a in records]), func)
        if not found:
            if (
                section == "pairs"
                and func == 1
                and self.top.defaults.gen_pairs == "yes"
            ):
                return [Interaction(section, atoms, func)]
            where = ", ".join(f"{a.name} ({a.type})" for a in records)
            raise TopologyError(
                f"No default [ {table} ] parameters (function {func}) for atoms "
                f"{where} of molecule {mol.name} at {line}."
            )

        b_types = [a.type_b or a.type for a in records]
        if any(a.type_b for a in records):
            b_found = self._index.find(table, names(b_types), func)
            if len(b_found) == len(found):
                return [
                    Interaction(
                        section,
                        atoms,
                        func,
                        {**a.params, **{f"{k}_b": v for k, v in b.params.items()}},
                    )
                    for a, b in zip(found, b_found)
                ]
        return [Interaction(section, atoms, func, dict(e.params)) for e in found]

    # -- system -----------------------------------------------------------

    def _system(self, line: Line, tokens: list[str]) -> None:
        self.top.system = (self.top.system + " " + line.text).strip()

    def _molecules(self, line: Line, tokens: list[str]) -> None:
        if len(tokens) != 2:
            raise TopologyError(f"[ molecules ] expects a name and a count at {line}.")
        try:
            mol = self.top.moleculetype(tokens[0])
        except KeyError:
            raise TopologyError(
                f"Unknown molecule type '{tokens[0]}' at {line}."
            ) from None
        self.top.molecules.append((mol, self._int(tokens[1], line, "molecule count")))


_HANDLERS = {
    "defaults": "_defaults",
    "atomtypes": "_atomtypes",
    "cmaptypes": "_cmaptypes",
    "moleculetype": "_moleculetype",
    "atoms": "_atoms",
    "exclusions": "_exclusions",
    "system": "_system",
    "molecules": "_molecules",
}
