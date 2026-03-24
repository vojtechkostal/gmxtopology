from dataclasses import dataclass
from typing import Callable, Dict, List, Literal, Optional, Tuple

from .interaction_specs import (
    ANGLE_RESTRAINTS,
    ANGLE_RESTRAINTS_Z,
    ANGLES,
    BONDS,
    CONSTRAINTS,
    DIHEDRAL_RESTRAINTS,
    DIHEDRALS,
    DISTANCE_RESTRAINTS,
    NONBOND_PARAMS,
    ORIENTATION_RESTRAINTS,
    PAIRS,
    PAIRS_NB,
    POSITION_RESTRAINTS,
    SETTLES,
    VIRTUAL_SITES1,
    VIRTUAL_SITES2,
    VIRTUAL_SITES3,
    VIRTUAL_SITES4,
    VIRTUAL_SITESN,
)
from .lookup import lookup_paramtype
from .schema import InteractionSpec
from .topology import (
    Angle,
    AngleRestraint,
    AngleRestraintZ,
    AngleType,
    Atom,
    AtomType,
    Bond,
    BondType,
    Constraint,
    ConstraintType,
    Defaults,
    Dihedral,
    DihedralRestraint,
    DihedralType,
    DistanceRestraint,
    Exclusion,
    MoleculeType,
    NonBondParam,
    OrientationRestraint,
    Pair,
    PairNB,
    PairType,
    PositionRestraint,
    Settle,
    System,
    Topology,
    VirtualSite1,
    VirtualSite2,
    VirtualSite3,
    VirtualSite4,
    VirtualSiteN,
    _find_atom,
    _find_atomtype,
    _find_molecule_type,
)

SectionScope = Literal["topology", "molecule"]


@dataclass(frozen=True)
class InteractionSection:
    scope: SectionScope
    spec: InteractionSpec
    record_type: type


def _parse_int(token: str, *, ctx: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise ValueError(f"Invalid integer value '{token}' in {ctx}.") from exc


def _resolve_atomtypes(tokens: List[str], top: Topology) -> list[AtomType]:
    return [_find_atomtype(top, token) for token in tokens]


def _resolve_atoms(tokens: List[str], mol: MoleculeType, ctx: str) -> list[Atom]:
    return [_find_atom(mol, token, ctx) for token in tokens]


def parse_molecule(
    parts: List[str],
    top: Topology,
) -> Tuple[MoleculeType, int]:
    return _find_molecule_type(top, parts[0]), int(parts[1])


def parse_interaction_type(
    parts: List[str],
    interaction_spec: InteractionSpec,
    section_cls: Callable,
    top: Topology,
    section: str = "",
    ifdef_state: Optional[str] = None,
) -> object:
    n_atom = interaction_spec.n_atoms
    atomtypes = _resolve_atomtypes(parts[:n_atom], top)

    func = _parse_int(parts[n_atom], ctx=f"section {section} in topology {top.source}")
    param_tokens = parts[n_atom + 1:]
    params = interaction_spec.parse(
        func,
        param_tokens,
        ctx=f"section {section} in topology {top.source}",
    )
    if not param_tokens:
        raise ValueError(
            f"Parameters must be provided for {section} in {top.source} "
            f"with atomtypes ({', '.join(atom.name for atom in atomtypes)})."
        )

    return section_cls(*atomtypes, func, params, ifdef_state=ifdef_state)


def parse_interaction(
    parts: List[str],
    interaction_spec: InteractionSpec,
    section_cls: object,
    top: Topology,
    mol: MoleculeType,
    section: str = "",
    ifdef_state: Optional[str] = None,
) -> list[Bond, Pair, Angle, Dihedral, Exclusion]:
    n_atoms = interaction_spec.n_atoms
    ctx = f"section {section} in molecule {mol.name}"
    atoms = _resolve_atoms(parts[:n_atoms], mol, ctx)
    func = _parse_int(parts[n_atoms], ctx=ctx)
    param_tokens = parts[n_atoms + 1:]

    if not param_tokens:
        if top.defaults.gen_pairs == "yes" and section in {"pairs", "pairs_nb"}:
            return [section_cls(*atoms, func, {})]
        paramtypes = lookup_paramtype(top, *atoms, func=func, section=section)
        return [section_cls(*atoms, func, item.params) for item in paramtypes]

    params = interaction_spec.parse(func, param_tokens, ctx=ctx)
    return [section_cls(*atoms, func, params, ifdef_state=ifdef_state)]


def apply_section_line(
    section: str,
    parts: List[str],
    top: Topology,
    active_molecule: Optional[MoleculeType],
    ifdef_state: str = "free",
) -> Optional[MoleculeType]:
    """Apply one parsed data line to the topology and return the active molecule."""

    section_entry = SECTION_REGISTRY.get(section)

    if section_entry is not None and section_entry.scope == "molecule" and (
        active_molecule is None
    ):
        raise ValueError(
            f"[ {section} ] section found before molecule type in {top.source}."
        )

    if section_entry is not None and section_entry.scope == "topology" and (
        active_molecule is not None
    ):
        raise ValueError(
            f"[ {section} ] section found inside molecule type "
            f"{active_molecule.name} in {top.source}. "
            f"Parameter types must be defined globally."
        )

    if section == "defaults":
        top.defaults = Defaults.from_line(parts, top, ifdef_state)
        return active_molecule

    if section == "atomtypes":
        top.atomtypes.append(AtomType.from_line(parts, top, ifdef_state))
        return active_molecule

    if section_entry is not None and section_entry.scope == "topology":
        interaction_type = parse_interaction_type(
            parts,
            section_entry.spec,
            section_entry.record_type,
            top,
            section,
            ifdef_state,
        )
        getattr(top, section).append(interaction_type)
        return active_molecule

    if section == "moleculetype":
        molecule = MoleculeType.from_line(parts, top, ifdef_state)
        top.moleculetypes.append(molecule)
        return molecule

    if section == "atoms":
        assert active_molecule is not None
        active_molecule.atoms.append(
            Atom.from_line(parts, top, active_molecule, ifdef_state)
        )
        return active_molecule

    if section_entry is not None and section_entry.scope == "molecule":
        assert active_molecule is not None
        interactions = parse_interaction(
            parts,
            section_entry.spec,
            section_entry.record_type,
            top,
            active_molecule,
            section,
            ifdef_state,
        )
        getattr(active_molecule, section).extend(interactions)
        return active_molecule

    if section == "exclusions":
        assert active_molecule is not None
        active_molecule.exclusions.append(
            Exclusion.from_line(parts, active_molecule, ifdef_state)
        )
        return active_molecule

    if section == "system":
        top.system = System.from_line(parts, top)
        return active_molecule

    if section == "molecules":
        molecule, count = parse_molecule(parts, top)
        top.molecules[molecule.name] = (molecule, count)
        return active_molecule

    return active_molecule


SECTION_REGISTRY: Dict[str, InteractionSection] = {
    "bondtypes": InteractionSection("topology", BONDS, BondType),
    "pairtypes": InteractionSection("topology", PAIRS, PairType),
    "angletypes": InteractionSection("topology", ANGLES, AngleType),
    "dihedraltypes": InteractionSection("topology", DIHEDRALS, DihedralType),
    "constrainttypes": InteractionSection(
        "topology",
        CONSTRAINTS,
        ConstraintType,
    ),
    "nonbond_params": InteractionSection(
        "topology",
        NONBOND_PARAMS,
        NonBondParam,
    ),
    "bonds": InteractionSection("molecule", BONDS, Bond),
    "pairs": InteractionSection("molecule", PAIRS, Pair),
    "pairs_nb": InteractionSection("molecule", PAIRS_NB, PairNB),
    "angles": InteractionSection("molecule", ANGLES, Angle),
    "dihedrals": InteractionSection("molecule", DIHEDRALS, Dihedral),
    "constraints": InteractionSection("molecule", CONSTRAINTS, Constraint),
    "settles": InteractionSection("molecule", SETTLES, Settle),
    "virtual_sites1": InteractionSection(
        "molecule",
        VIRTUAL_SITES1,
        VirtualSite1,
    ),
    "virtual_sites2": InteractionSection(
        "molecule",
        VIRTUAL_SITES2,
        VirtualSite2,
    ),
    "virtual_sites3": InteractionSection(
        "molecule",
        VIRTUAL_SITES3,
        VirtualSite3,
    ),
    "virtual_sites4": InteractionSection(
        "molecule",
        VIRTUAL_SITES4,
        VirtualSite4,
    ),
    "virtual_sitesn": InteractionSection(
        "molecule",
        VIRTUAL_SITESN,
        VirtualSiteN,
    ),
    "dummies1": InteractionSection("molecule", VIRTUAL_SITES1, VirtualSite1),
    "dummies2": InteractionSection("molecule", VIRTUAL_SITES2, VirtualSite2),
    "dummies3": InteractionSection("molecule", VIRTUAL_SITES3, VirtualSite3),
    "dummies4": InteractionSection("molecule", VIRTUAL_SITES4, VirtualSite4),
    "dummiesn": InteractionSection("molecule", VIRTUAL_SITESN, VirtualSiteN),
    "position_restraints": InteractionSection(
        "molecule",
        POSITION_RESTRAINTS,
        PositionRestraint,
    ),
    "distance_restraints": InteractionSection(
        "molecule",
        DISTANCE_RESTRAINTS,
        DistanceRestraint,
    ),
    "dihedral_restraints": InteractionSection(
        "molecule",
        DIHEDRAL_RESTRAINTS,
        DihedralRestraint,
    ),
    "orientation_restraints": InteractionSection(
        "molecule",
        ORIENTATION_RESTRAINTS,
        OrientationRestraint,
    ),
    "angle_restraints": InteractionSection(
        "molecule",
        ANGLE_RESTRAINTS,
        AngleRestraint,
    ),
    "angle_restraints_z": InteractionSection(
        "molecule",
        ANGLE_RESTRAINTS_Z,
        AngleRestraintZ,
    ),
}

MOLECULE_SECTIONS: tuple[str, ...] = tuple(
    name for name, entry in SECTION_REGISTRY.items() if entry.scope == "molecule"
)
