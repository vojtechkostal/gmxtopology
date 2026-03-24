from __future__ import annotations

import types
from pathlib import Path
from dataclasses import dataclass, field
from typing import (
    Any,
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
    Set
)


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

    def _field_type(self, annotation: Any) -> Optional[type]:
        if annotation in {float, int, str}:
            return annotation

        origin = get_origin(annotation)
        if origin not in (types.UnionType, Union):
            return None

        inner_types = [
            arg for arg in get_args(annotation)
            if isinstance(arg, type) and arg is not type(None)
        ]
        if len(inner_types) == 1 and inner_types[0] in {float, int, str}:
            return inner_types[0]
        return None

    def _render_value(
        self,
        field_name: str,
        annotation: Any,
        value: Any,
    ) -> Optional[tuple[Any, type]]:
        scalar_type = self._field_type(annotation)
        if scalar_type is not None:
            return value, scalar_type

        if annotation in {AtomType, MoleculeType}:
            return value.name, str

        if annotation == Atom:
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

    def _fmt(self, name, type_: type) -> int:
        """Determine format string for a given parameter based on its name and type."""
        if type_ == float:
            width = self.W_FLOAT
            if name == "charge":
                decimals = self.D_CHARGE
            elif name in ("sigma", "epsilon"):
                decimals = self.D_LJ
            else:
                decimals = self.D

            fmt = f'>{width}.{decimals}f'
        else:
            width = self.W
            fmt = f'>{width}'

        return fmt

    def _render(self, *, values: bool) -> str:
        """Render either the legend (parameter names) or values."""

        chunks: list[str] = []
        for i, (key, (val, typ)) in enumerate(self.args.items()):
            if i == 0:
                if values:
                    chunks.append(f"{val:>{self.W_FIRST}}")
                else:
                    chunks.append(f"; {key:>{self.W_FIRST - 2}}")
                continue

            fmt = self._fmt(key, typ)
            if values:
                try:
                    chunks.append(f"{val:{fmt}}")
                except ValueError:
                    # Fallback for non-numeric values
                    chunks.append(f"{str(val):>{self.W}}")
            else:
                fmt = fmt.split('.')[0]  # remove decimal part for legend
                chunks.append(f"{key:{fmt}}")

        return "".join(chunks)

    @property
    def ifdef(self) -> Optional[str]:
        """Return the ifdef directive string for this section."""
        if self.ifdef_state:
            if "ifdef" in self.ifdef_state:
                return f"#{self.ifdef_state}"
            elif "else" in self.ifdef_state:
                return "#else"
        return self.ifdef_state

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
    ifdef_state: Optional[str] = "free"

    SECTION = 'defaults'


@dataclass(slots=True)
class Define(UpdatableSection):
    directive: str
    argument: Optional[str | int | float]

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"argument"})

    def __hash__(self):
        return hash(self.directive)

    def __str__(self) -> str:
        if self.argument is not None:
            return f"#define {self.directive} {self.argument}"
        return f"#define {self.directive}"


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
    ifdef_state: Optional[str] = "free"

    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"sigma", "epsilon"})
    SECTION = 'atomtypes'

    def __hash__(self):
        return hash(self.name)


@dataclass(slots=True)
class BondType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'bondtypes'


@dataclass(slots=True)
class PairType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairtypes'


@dataclass(slots=True)
class AngleType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angletypes'


@dataclass(slots=True)
class DihedralType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    ak: AtomType
    al: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedraltypes'


@dataclass(slots=True)
class ConstraintType(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'constrainttypes'


@dataclass(slots=True)
class Constraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'constraints'


@dataclass(slots=True)
class NonBondParam(UpdatableSection):
    ai: AtomType
    aj: AtomType
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

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
    mass: float
    ifdef_state: Optional[str] = "free"

    SECTION = 'atoms'
    MODIFIABLE: ClassVar[FrozenSet[str]] = frozenset({"charge"})


@dataclass(slots=True)
class Bond(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'bonds'


@dataclass(slots=True)
class Pair(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairs'


@dataclass(slots=True)
class PairNB(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'pairs_nb'


@dataclass(slots=True)
class Angle(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angles'


@dataclass(slots=True)
class Dihedral(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedrals'

    def __repr__(self):
        return (
            "Dihedral between atoms "
            f"{self.ai.name}, {self.aj.name}, {self.ak.name}, {self.al.name}"
        )


@dataclass(slots=True)
class Exclusion(RenderableSection):
    excluded: List[Atom]  # variable length per line
    ifdef_state: Optional[str] = "free"

    SECTION = 'exclusions'


@dataclass(slots=True)
class Settle(UpdatableSection):
    ai: Atom
    func: int
    params: Dict[str, float] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'settles'


@dataclass(slots=True)
class VirtualSite1(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites1'


@dataclass(slots=True)
class VirtualSite2(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites2'


@dataclass(slots=True)
class VirtualSite3(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

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
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sites4'


@dataclass(slots=True)
class VirtualSiteN(UpdatableSection):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'virtual_sitesn'


@dataclass(slots=True)
class PositionRestraint(UpdatableSection):
    ai: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'position_restraints'


@dataclass(slots=True)
class DistanceRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'distance_restraints'


@dataclass(slots=True)
class DihedralRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'dihedral_restraints'


@dataclass(slots=True)
class OrientationRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'orientation_restraints'


@dataclass(slots=True)
class AngleRestraint(UpdatableSection):
    ai: Atom
    aj: Atom
    ak: Atom
    al: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angle_restraints'


@dataclass(slots=True)
class AngleRestraintZ(UpdatableSection):
    ai: Atom
    aj: Atom
    func: int
    params: Mapping[str, float | int | str] = field(default_factory=dict)
    ifdef_state: Optional[str] = "free"

    SECTION = 'angle_restraints_z'


@dataclass(slots=True)
class System(RenderableSection):
    description: str = "system"

    SECTION = 'system'


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
    ifdef_state: Optional[str] = "free"

    SECTION = 'moleculetype'

    def get_atom_by_idx(self, idx: int) -> Optional[Atom]:
        for atom in self.atoms:
            if atom.nr == idx:
                return atom
        return None

    @property
    def _vsite_sections(self):
        return {
            attr: getattr(self, attr)
            for attr in self.__dataclass_fields__
            if attr.startswith("virtual_sites")
        }

    @property
    def _connection_sections(self):
        return {
            "bonds", "pairs", "angles", "dihedrals", "constraints",
            "position_restraints", "distance_restraints",
            "dihedral_restraints", "orientation_restraints",
            "angle_restraints", "angle_restraints_z"
        }

    def remove_vsites(self) -> None:
        from .parser import remove_vsites
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
    moleculetypes: List[MoleculeType] = field(default_factory=list)

    system: Optional[System] = None
    molecules: Dict[str, Tuple[MoleculeType, int]] = field(default_factory=dict)

    defines: Set[Define] = field(default_factory=set)

    @property
    def residues(self) -> List[MoleculeType]:
        """List of all residues in the topology, expanded by their counts."""
        return [mol for mol, count in self.molecules.values() for _ in range(count)]

    @property
    def atoms(self) -> List[Atom]:
        """List of all atoms in the topology."""
        return [atom for mol in self.residues for atom in mol.atoms]

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
