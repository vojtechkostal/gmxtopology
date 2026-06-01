from __future__ import annotations

import types
from pathlib import Path
from dataclasses import dataclass, field
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    ClassVar,
    FrozenSet,
    Mapping,
    get_type_hints,
    get_origin,
    get_args,
    Union,
    Tuple,
)

PreprocessorState = tuple[str, ...]


def _parse_int(token: str, *, ctx: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise ValueError(f"Invalid integer value '{token}' in {ctx}.") from exc


def _parse_float(token: str, *, ctx: str) -> float:
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(f"Invalid float value '{token}' in {ctx}.") from exc


def _parse_define_or_float(
    token: str,
    top: Topology,
    *,
    ctx: str,
) -> float | str:
    try:
        return float(token)
    except ValueError:
        define_names = {define.directive for define in top.defines}
        if token.strip("-") in define_names:
            return token
        raise ValueError(f"Invalid numeric value '{token}' in {ctx}.")


def _ensure_unique(
    items: list[Any],
    predicate: Callable[[Any], bool],
    message: str,
) -> None:
    if any(predicate(item) for item in items):
        raise ValueError(message)


def _find_atomtype(top: Topology, name: str) -> AtomType:
    if name == "X":
        return AtomType(
            name="X",
            bonded_type=None,
            atnum=0,
            mass=0.0,
            charge=0.0,
            ptype="A",
            sigma=0.0,
            epsilon=0.0,
        )

    cached = getattr(top, "_atomtype_index", None)
    if cached is None or cached[0] != len(top.atomtypes):
        atomtypes = {}
        for atomtype in top.atomtypes:
            atomtypes.setdefault(atomtype.name, atomtype)
        cached = len(top.atomtypes), atomtypes
        top._atomtype_index = cached

    try:
        return cached[1][name]
    except KeyError as exc:
        raise ValueError(f"Atomtype {name} not found in {top.source}.") from exc


def _find_molecule_type(top: Topology, name: str) -> MoleculeType:
    for moleculetype in top.moleculetypes:
        if moleculetype.name.casefold() == name.casefold():
            return moleculetype

    raise ValueError(
        f"Molecule type '{name}' not found in topology {top.source}."
    )


def _find_atom(mol: MoleculeType, token: str, ctx: str) -> Atom:
    if not token.isdigit():
        raise ValueError(f"Invalid atom token '{token}' in {ctx}.")

    atom = mol.get_atom_by_idx(int(token))
    if atom is None:
        raise ValueError(f"Atom index '{token}' not found in {ctx}.")

    return atom


class RenderableSection:
    """Shared rendering helpers for topology records."""

    SECTION = ClassVar[str]
    W = 10
    W_FIRST = 8
    W_FLOAT = 15
    D = 4
    D_CHARGE = 10
    D_LJ = 6

    @property
    def header(self) -> str:
        """Section header string."""
        return f"[ {self.SECTION} ]"

    @staticmethod
    def _unwrap_optional(annotation: Any) -> Any:
        origin = get_origin(annotation)
        if origin in (types.UnionType, Union):
            inner_types = [arg for arg in get_args(annotation) if arg is not type(None)]
            if len(inner_types) == 1:
                return inner_types[0]
        return annotation

    def _render_value(
        self,
        field_name: str,
        annotation: Any,
        value: Any,
    ) -> Optional[tuple[Any, type]]:
        field_type = self._unwrap_optional(annotation)
        if field_type in {float, int, str}:
            return value, field_type

        if field_type in {AtomType, MoleculeType}:
            return value.name, str

        if field_type == Atom:
            return value.nr, int

        if field_name == "excluded":
            return " ".join(str(atom.nr) for atom in value), str

        return None

    @property
    def args(self) -> Dict[str, Any]:
        """Extract fields and their display values for rendering."""
        hints = get_type_hints(type(self))
        out: Dict[str, Tuple[Any, type]] = {}

        for field_name, annotation in hints.items():
            if field_name in {"MODIFIABLE", "SECTION", "ifdef_state"}:
                continue

            value = getattr(self, field_name)
            if value is None:
                continue

            if field_name == "params":
                out.update({
                    param_name: (param_value, type(param_value))
                    for param_name, param_value in value.items()
                })
                continue

            rendered = self._render_value(field_name, annotation, value)
            if rendered is not None:
                out[field_name] = rendered

        return out

    def _fmt(self, name: str, type_: type) -> str:
        """Determine format string for a given parameter based on its name and type."""
        if type_ == float:
            if name == "charge":
                decimals = self.D_CHARGE
            elif name in ("sigma", "epsilon"):
                decimals = self.D_LJ
            else:
                decimals = self.D
            return f'>{self.W_FLOAT}.{decimals}f'
        return f'>{self.W}'

    def _render(self, *, values: bool) -> str:
        """Render either the legend (parameter names) or values."""

        chunks: list[str] = []
        for i, (name, (value, type_)) in enumerate(self.args.items()):
            if i == 0:
                if values:
                    chunks.append(f"{value:>{self.W_FIRST}}")
                else:
                    chunks.append(f"; {name:>{self.W_FIRST - 2}}")
                continue

            fmt = self._fmt(name, type_)
            if values:
                try:
                    chunks.append(f"{value:{fmt}}")
                except ValueError:
                    # Fallback for non-numeric values
                    chunks.append(f" {str(value):>{self.W}}")
            else:
                fmt = fmt.split('.')[0]  # remove decimal part for legend
                chunks.append(f"{name:{fmt}}")

        return "".join(chunks)

    @property
    def legend(self) -> str:
        """Return the legend string (parameter names)."""
        return self._render(values=False)

    def __str__(self) -> str:
        """Return the values string (parameter values)."""
        return self._render(values=True)


class UpdatableSection(RenderableSection):
    """Renderable section that allows restricted in-place updates."""

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset()

    def update(self, **kwargs: Any) -> None:
        params = getattr(self, "params", None)
        if isinstance(params, dict):
            unknown = kwargs.keys() - params.keys()
            if unknown:
                raise KeyError(
                    f"Cannot add new parameter(s) {sorted(unknown)}; "
                    f"allowed keys are {sorted(params.keys())}"
                )
            params.update(kwargs)
            return

        unknown = kwargs.keys() - type(self).MODIFIABLE
        if unknown:
            raise KeyError(
                f"Parameter(s) {sorted(unknown)} are not modifiable "
                f"for {type(self).__name__}; allowed are "
                f"{sorted(type(self).MODIFIABLE)}"
            )
        for key, value in kwargs.items():
            setattr(self, key, value)


@dataclass(slots=True)
class Defaults(RenderableSection):
    nbfunc: int
    comb_rule: int
    gen_pairs: str
    fudgeLJ: float
    fudgeQQ: float
    n: Optional[int] = None
    ifdef_state: PreprocessorState = ()

    SECTION = 'defaults'

    @classmethod
    def from_line(
        cls,
        parts: list[str],
        top: Topology,
        ifdef_state: PreprocessorState = (),
    ) -> Defaults:
        ctx = f"[ defaults ] in {top.source}"
        if top.defaults is not None:
            raise ValueError("Multiple [ defaults ] sections found.")
        if len(parts) < 2 or len(parts) > 6:
            raise ValueError(f"{ctx} expects 2 to 6 values, got {len(parts)}.")

        defaults = cls(
            nbfunc=_parse_int(parts[0], ctx=ctx),
            comb_rule=_parse_int(parts[1], ctx=ctx),
            gen_pairs=parts[2] if len(parts) >= 3 else "no",
            fudgeLJ=_parse_float(parts[3], ctx=ctx) if len(parts) >= 4 else 1.0,
            fudgeQQ=_parse_float(parts[4], ctx=ctx) if len(parts) >= 5 else 1.0,
            n=_parse_int(parts[5], ctx=ctx) if len(parts) == 6 else None,
            ifdef_state=ifdef_state,
        )

        if defaults.nbfunc not in (1, 2):
            raise ValueError(
                f"Invalid nbfunc in {ctx}: expected 1 or 2, "
                f"got {defaults.nbfunc}"
            )

        if defaults.comb_rule not in (1, 2, 3):
            raise ValueError(
                f"Invalid comb_rule in {ctx}: expected 1, 2, or 3, "
                f"got {defaults.comb_rule}"
            )

        if defaults.gen_pairs not in {"yes", "no"}:
            raise ValueError(
                f"Invalid gen_pairs in {ctx}: expected 'yes' or 'no', "
                f"got '{defaults.gen_pairs}'"
            )

        return defaults


@dataclass(slots=True)
class Define(UpdatableSection):
    directive: str
    argument: Optional[str | int | float]
    ifdef_state: PreprocessorState = ()

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"argument"})

    def __str__(self) -> str:
        if self.argument is not None:
            return f"#define {self.directive} {self.argument}"
        return f"#define {self.directive}"


@dataclass(slots=True)
class RawSection:
    """Topology section preserved without modeling its contents."""

    name: str
    lines: List[str] = field(default_factory=list)
    ifdef_state: PreprocessorState = ()

    @property
    def header(self) -> str:
        return f"[ {self.name} ]"


@dataclass(slots=True)
class AtomType(UpdatableSection):
    name: str
    bonded_type: Optional[str]
    atnum: Optional[int]
    mass: float
    charge: float
    ptype: str
    sigma: float
    epsilon: float
    ifdef_state: PreprocessorState = ()

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"sigma", "epsilon"})
    SECTION = 'atomtypes'

    def __hash__(self):
        return hash(self.name)

    @classmethod
    def from_line(
        cls,
        parts: list[str],
        top: Topology,
        ifdef_state: PreprocessorState = (),
    ) -> AtomType:
        ctx = f"[ atomtypes ] in {top.source}"
        if len(parts) not in {6, 7, 8}:
            raise ValueError(f"{ctx} expects 6, 7, or 8 values, got {len(parts)}.")

        bonded_type: str | None = None
        atnum: int | None = None
        mass_idx = 1

        if len(parts) == 7 and parts[1].isdigit():
            atnum = _parse_int(parts[1], ctx=ctx)
            mass_idx = 2
        elif len(parts) == 7:
            bonded_type = parts[1]
            mass_idx = 2
        elif len(parts) == 8:
            bonded_type = parts[1]
            atnum = _parse_int(parts[2], ctx=ctx)
            mass_idx = 3

        atomtype = cls(
            name=parts[0],
            bonded_type=bonded_type,
            atnum=atnum,
            mass=_parse_define_or_float(parts[mass_idx], top, ctx=ctx),
            charge=_parse_define_or_float(parts[mass_idx + 1], top, ctx=ctx),
            ptype=parts[mass_idx + 2],
            sigma=_parse_define_or_float(parts[mass_idx + 3], top, ctx=ctx),
            epsilon=_parse_define_or_float(parts[mass_idx + 4], top, ctx=ctx),
            ifdef_state=ifdef_state,
        )

        if atomtype.ptype not in {"A", "S", "V", "D"}:
            raise ValueError(
                f"Invalid ptype of atomtype {atomtype.name} in {ctx}: "
                f"expected 'A', 'S', 'V', or 'D', got '{atomtype.ptype}'"
            )

        _ensure_unique(
            top.atomtypes,
            lambda item: item.name == atomtype.name
            and item.ifdef_state == atomtype.ifdef_state,
            f"Duplicate atomtype '{atomtype.name}' found in topology "
            f"{top.source}.",
        )
        return atomtype


@dataclass(slots=True)
class BondType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'bondtypes'


@dataclass(slots=True)
class PairType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'pairtypes'


@dataclass(slots=True)
class AngleType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'angletypes'


@dataclass(slots=True)
class DihedralType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    al: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'dihedraltypes'


@dataclass(slots=True)
class ConstraintType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'constrainttypes'


@dataclass(slots=True)
class Constraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'constraints'


@dataclass(slots=True)
class NonBondParam(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'nonbond_params'

    def __hash__(self):
        return hash((self.ai.name, self.aj.name))


@dataclass(slots=True)
class Atom(UpdatableSection):
    nr: int
    type: AtomType
    resnr: int
    residue: MoleculeType
    name: str
    cgnr: int
    charge: float
    mass: Optional[float]
    type_b: Optional[AtomType] = None
    charge_b: Optional[float] = None
    mass_b: Optional[float] = None
    ifdef_state: PreprocessorState = ()

    SECTION = 'atoms'
    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"charge", "charge_b"})

    @classmethod
    def from_line(
        cls,
        parts: list[str],
        top: Topology,
        mol: MoleculeType,
        ifdef_state: PreprocessorState = (),
    ) -> Atom:
        ctx = f"[ atoms ] in {top.source}"
        if len(parts) not in {7, 8, 10, 11}:
            raise NotImplementedError(
                f"{ctx} expects 7 values, an optional mass, and optionally "
                "the three topology-B values; "
                f"got {len(parts)}."
            )

        nr = _parse_int(parts[0], ctx=ctx)
        if mol.get_atom_by_idx(nr) is not None and not ifdef_state:
            raise ValueError(
                f"Duplicate atom index '{nr}' in molecule '{mol.name}' "
                f"in topology {top.source}"
            )

        has_mass = len(parts) in {8, 11}
        has_b_state = len(parts) in {10, 11}
        b_state_idx = 8 if has_mass else 7

        return cls(
            nr=nr,
            type=_find_atomtype(top, parts[1]),
            resnr=_parse_int(parts[2], ctx=ctx),
            residue=mol,
            name=parts[4],
            cgnr=_parse_int(parts[5], ctx=ctx),
            charge=_parse_define_or_float(parts[6], top, ctx=ctx),
            mass=(
                _parse_define_or_float(parts[7], top, ctx=ctx)
                if has_mass
                else None
            ),
            type_b=(
                _find_atomtype(top, parts[b_state_idx])
                if has_b_state
                else None
            ),
            charge_b=(
                _parse_define_or_float(parts[b_state_idx + 1], top, ctx=ctx)
                if has_b_state
                else None
            ),
            mass_b=(
                _parse_define_or_float(parts[b_state_idx + 2], top, ctx=ctx)
                if has_b_state
                else None
            ),
            ifdef_state=ifdef_state,
        )


@dataclass(slots=True)
class Bond(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'bonds'


@dataclass(slots=True)
class Pair(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'pairs'


@dataclass(slots=True)
class PairNB(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'pairs_nb'


@dataclass(slots=True)
class Angle(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'angles'


@dataclass(slots=True)
class Dihedral(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'dihedrals'

    def __repr__(self):
        return (
            "Dihedral between atoms "
            f"{self.ai.name}, {self.aj.name}, {self.ak.name}, {self.al.name}"
        )


@dataclass(slots=True)
class Exclusion(RenderableSection):
    excluded: List[Atom]  # variable length per line
    ifdef_state: PreprocessorState = ()

    SECTION = 'exclusions'

    @classmethod
    def from_line(
        cls,
        parts: list[str],
        mol: MoleculeType,
        ifdef_state: PreprocessorState = (),
    ) -> Exclusion:
        excluded = [_find_atom(mol, token, "exclusions") for token in parts]
        return cls(excluded, ifdef_state=ifdef_state)


@dataclass(slots=True)
class Settle(UpdatableSection):
    ai: Atom
    func: int
    params: Dict[str, float] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'settles'


@dataclass(slots=True)
class VirtualSite1(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'virtual_sites1'


@dataclass(slots=True)
class VirtualSite2(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'virtual_sites2'


@dataclass(slots=True)
class VirtualSite3(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'virtual_sites3'


@dataclass(slots=True)
class VirtualSite4(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    am: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'virtual_sites4'


@dataclass(slots=True)
class VirtualSiteN(UpdatableSection):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'virtual_sitesn'


@dataclass(slots=True)
class PositionRestraint(UpdatableSection):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'position_restraints'


@dataclass(slots=True)
class DistanceRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'distance_restraints'


@dataclass(slots=True)
class DihedralRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'dihedral_restraints'


@dataclass(slots=True)
class OrientationRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'orientation_restraints'


@dataclass(slots=True)
class AngleRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'angle_restraints'


@dataclass(slots=True)
class AngleRestraintZ(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: PreprocessorState = ()

    SECTION = 'angle_restraints_z'


@dataclass(slots=True)
class System(RenderableSection):
    description: str = "system"

    SECTION = 'system'

    @classmethod
    def from_line(cls, parts: list[str], top: Topology) -> System:
        if top.system is not None:
            raise ValueError("Multiple [ system ] sections found.")
        description = " ".join(parts) if parts else "system"
        return cls(description)


@dataclass(slots=True)
class MoleculeType(RenderableSection):
    name: str
    nrexcl: int
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Bond] = field(default_factory=list)
    pairs: List[Pair] = field(default_factory=list)
    pairs_nb: List[PairNB] = field(default_factory=list)
    angles: List[Angle] = field(default_factory=list)
    dihedrals: List[Dihedral] = field(default_factory=list)
    constraints: List[Constraint] = field(default_factory=list)
    exclusions: List[Exclusion] = field(default_factory=list)
    settles: List[Settle] = field(default_factory=list)
    virtual_sites1: List[VirtualSite1] = field(default_factory=list)
    virtual_sites2: List[VirtualSite2] = field(default_factory=list)
    virtual_sites3: List[VirtualSite3] = field(default_factory=list)
    virtual_sites4: List[VirtualSite4] = field(default_factory=list)
    virtual_sitesn: List[VirtualSiteN] = field(default_factory=list)
    position_restraints: List[PositionRestraint] = field(default_factory=list)
    distance_restraints: List[DistanceRestraint] = field(default_factory=list)
    dihedral_restraints: List[DihedralRestraint] = field(default_factory=list)
    orientation_restraints: List[OrientationRestraint] = field(default_factory=list)
    angle_restraints: List[AngleRestraint] = field(default_factory=list)
    angle_restraints_z: List[AngleRestraintZ] = field(default_factory=list)
    raw_sections: List[RawSection] = field(default_factory=list)
    ifdef_state: PreprocessorState = ()

    SECTION = 'moleculetype'
    VSITE_SECTIONS = (
        "virtual_sites1",
        "virtual_sites2",
        "virtual_sites3",
        "virtual_sites4",
        "virtual_sitesn",
    )
    CONNECTION_SECTIONS = (
        "bonds",
        "pairs",
        "angles",
        "dihedrals",
        "constraints",
        "position_restraints",
        "distance_restraints",
        "dihedral_restraints",
        "orientation_restraints",
        "angle_restraints",
        "angle_restraints_z",
    )

    @classmethod
    def from_line(
        cls,
        parts: list[str],
        top: Topology,
        ifdef_state: PreprocessorState = (),
    ) -> MoleculeType:
        molecule = cls(
            name=parts[0],
            nrexcl=_parse_int(parts[1], ctx=f"[ moleculetype ] in {top.source}"),
            ifdef_state=ifdef_state,
        )
        _ensure_unique(
            top.moleculetypes,
            lambda item: item.name.casefold() == molecule.name.casefold()
            and item.ifdef_state == molecule.ifdef_state,
            f"Duplicate moleculetype '{molecule.name}' found in topology "
            f"{top.source}.",
        )
        return molecule

    def get_atom_by_idx(self, idx: int) -> Optional[Atom]:
        for atom in self.atoms:
            if atom.nr == idx:
                return atom
        return None

    def remove_vsites(self) -> None:
        from .lookup import remove_vsites
        remove_vsites(self)

    def __repr__(self) -> str:
        return f"<MoleculeType name={self.name} nrexcl={self.nrexcl}>"


@dataclass(slots=True)
class Topology:
    """Representation of a molecular topology file."""
    source: Path
    defaults: Optional[Defaults] = None
    atomtypes: List[AtomType] = field(default_factory=list)
    nonbond_params: List[NonBondParam] = field(default_factory=list)
    bondtypes: List[BondType] = field(default_factory=list)
    pairtypes: List[PairType] = field(default_factory=list)
    angletypes: List[AngleType] = field(default_factory=list)
    dihedraltypes: List[DihedralType] = field(default_factory=list)
    constrainttypes: List[ConstraintType] = field(default_factory=list)
    raw_sections: List[RawSection] = field(default_factory=list)
    moleculetypes: List[MoleculeType] = field(default_factory=list)

    system: Optional[System] = None
    molecules: Dict[str, Tuple[MoleculeType, int]] = field(default_factory=dict)

    defines: List[Define] = field(default_factory=list)

    _atomtype_index: Optional[Tuple[int, Dict[str, AtomType]]] = field(
        default=None,
        init=False,
        repr=False,
    )
    _paramtype_indexes: Dict[str, Any] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )

    @property
    def residues(self) -> List[MoleculeType]:
        """List of all residues in the topology, expanded by their counts."""
        return [mol for mol, count in self.molecules.values() for _ in range(count)]

    @property
    def atoms(self) -> List[Atom]:
        """List of all atoms in the topology."""
        return [
            atom
            for mol, count in self.molecules.values()
            for _ in range(count)
            for atom in mol.atoms
        ]

    def __post_init__(self) -> None:
        self.source = Path(self.source).resolve()
        from .io import read_topology
        read_topology(self.source, self)

    def write(self, fn: str | Path, **kwargs) -> None:
        from .io import write_topology
        write_topology(self, Path(fn).resolve(), **kwargs)

    def __repr__(self) -> str:
        system_name = self.system.description if self.system is not None else "unset"
        molecule_count = sum(count for _, count in self.molecules.values())
        return (
            f"<Topology system='{system_name}' "
            f"moleculetypes={len(self.moleculetypes)} "
            f"molecules={molecule_count} "
            f"atoms={len(self.atoms)}>"
        )
