from typing import List, Dict, Tuple, Any, Callable, Optional

from ..topology import (
    Defaults,
    AtomType, BondType, PairType, AngleType, DihedralType, ConstraintType, NonBondParam,
    MoleculeType,
    Atom, Bond, Pair, PairNB, Angle, Dihedral, Constraint,
    Exclusion, Settle, System,
    VirtualSite1, VirtualSite2, VirtualSite3, VirtualSite4, VirtualSiteN,
    PositionRestraint, DistanceRestraint, DihedralRestraint,
    OrientationRestraint, AngleRestraint, AngleRestraintZ,
    Topology
)

from .helpers import lookup_paramtype
from ..schema import InteractionSpec
from ..interaction_specs import (
    BONDS, PAIRS, PAIRS_NB, ANGLES, DIHEDRALS,
    CONSTRAINTS, SETTLES, NONBOND_PARAMS,
    VIRTUAL_SITES1, VIRTUAL_SITES2, VIRTUAL_SITES3, VIRTUAL_SITES4, VIRTUAL_SITESN,
    POSITION_RESTRAINTS, DISTANCE_RESTRAINTS, DIHEDRAL_RESTRAINTS,
    ORIENTATION_RESTRAINTS, ANGLE_RESTRAINTS, ANGLE_RESTRAINTS_Z
)


SectionSpec = tuple[InteractionSpec, type]


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


def _parse_define_or_float(token: str, top: Topology, *, ctx: str) -> float | str:
    try:
        return float(token)
    except ValueError:
        define_names = {define.directive for define in top.defines}
        if token.strip("-") in define_names:
            return token
        raise ValueError(f"Invalid numeric value '{token}' in {ctx}.")


def _parse_defaults_values(parts: List[str], top: Topology) -> dict[str, Any]:
    if len(parts) < 2 or len(parts) > 6:
        raise ValueError(
            f"[ defaults ] in {top.source} expects 2 to 6 values, got {len(parts)}."
        )

    values = {
        "nbfunc": _parse_int(parts[0], ctx=f"[ defaults ] in {top.source}"),
        "comb_rule": _parse_int(parts[1], ctx=f"[ defaults ] in {top.source}"),
        "gen_pairs": parts[2] if len(parts) >= 3 else "no",
        "fudgeLJ": _parse_float(parts[3], ctx=f"[ defaults ] in {top.source}")
        if len(parts) >= 4 else 1.0,
        "fudgeQQ": _parse_float(parts[4], ctx=f"[ defaults ] in {top.source}")
        if len(parts) >= 5 else 1.0,
        "n": _parse_int(parts[5], ctx=f"[ defaults ] in {top.source}")
        if len(parts) == 6 else None,
    }
    return values


def _parse_atomtype_values(parts: List[str], top: Topology) -> dict[str, Any]:
    ctx = f"[ atomtypes ] in {top.source}"
    if len(parts) not in {6, 7, 8}:
        raise ValueError(f"{ctx} expects 6, 7, or 8 values, got {len(parts)}.")

    values: dict[str, Any] = {"name": parts[0], "bonded_type": None, "atnum": None}
    if len(parts) == 6:
        mass_idx = 1
    elif len(parts) == 7 and parts[1].isdigit():
        values["atnum"] = _parse_int(parts[1], ctx=ctx)
        mass_idx = 2
    elif len(parts) == 7:
        values["bonded_type"] = parts[1]
        mass_idx = 2
    else:
        values["bonded_type"] = parts[1]
        values["atnum"] = _parse_int(parts[2], ctx=ctx)
        mass_idx = 3

    values["mass"] = _parse_define_or_float(parts[mass_idx], top, ctx=ctx)
    values["charge"] = _parse_define_or_float(parts[mass_idx + 1], top, ctx=ctx)
    values["ptype"] = parts[mass_idx + 2]
    values["sigma"] = _parse_define_or_float(parts[mass_idx + 3], top, ctx=ctx)
    values["epsilon"] = _parse_define_or_float(parts[mass_idx + 4], top, ctx=ctx)
    return values


def _parse_atom_values(parts: List[str], top: Topology) -> dict[str, Any]:
    ctx = f"[ atoms ] in {top.source}"
    if len(parts) != 8:
        raise NotImplementedError(
            f"{ctx} currently supports exactly 8 values per atom line; "
            f"got {len(parts)}."
        )

    return {
        "nr": _parse_int(parts[0], ctx=ctx),
        "type": parts[1],
        "resnr": _parse_int(parts[2], ctx=ctx),
        "residue": parts[3],
        "name": parts[4],
        "cgnr": _parse_int(parts[5], ctx=ctx),
        "charge": _parse_define_or_float(parts[6], top, ctx=ctx),
        "mass": _parse_define_or_float(parts[7], top, ctx=ctx),
    }


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

    for atomtype in top.atomtypes:
        if atomtype.name == name:
            return atomtype

    raise ValueError(f"Atomtype {name} not found in {top.source}.")


def _find_molecule_type(top: Topology, name: str) -> MoleculeType:
    for moleculetype in top.moleculetypes:
        if moleculetype.name == name:
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


def _ensure_unique(
    items: list[Any],
    predicate: Callable[[Any], bool],
    message: str,
) -> None:
    if any(predicate(item) for item in items):
        raise ValueError(message)


def _resolve_atomtypes(tokens: List[str], top: Topology) -> list[AtomType]:
    return [_find_atomtype(top, token) for token in tokens]


def _resolve_atoms(tokens: List[str], mol: MoleculeType, ctx: str) -> list[Atom]:
    return [_find_atom(mol, token, ctx) for token in tokens]


def parse_defaults(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None
) -> Defaults:

    # check for defaults duplicates
    if top.defaults is not None:
        raise ValueError("Multiple [ defaults ] sections found.")

    inputs = _parse_defaults_values(parts, top)
    inputs["ifdef_state"] = ifdef_state

    if inputs["nbfunc"] not in (1, 2):
        raise ValueError(
            f"Invalid nbfunc in [ defaults ] in {top.source}: "
            f"expected 1 or 2, got {inputs['nbfunc']}"
        )

    if inputs["comb_rule"] not in (1, 2, 3):
        raise ValueError(
            f"Invalid comb_rule in [ defaults ] in {top.source}: "
            f"expected 1, 2, or 3, got {inputs['comb_rule']}"
        )

    if inputs["gen_pairs"] not in ('yes', 'no'):
        raise ValueError(
            f"Invalid gen_pairs in [ defaults ] in {top.source}: "
            f"expected 'yes' or 'no', got '{inputs['gen_pairs']}'"
        )

    return Defaults(**inputs)


def parse_atomtype(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None,
) -> AtomType:

    inputs = _parse_atomtype_values(parts, top)
    inputs["ifdef_state"] = ifdef_state

    if inputs['ptype'] not in ("A", "S", "V", "D"):
        raise ValueError(
            f"Invalid ptype of atomtype {inputs['name']} "
            f"in [ atomtypes ] in {top.source}: "
            f"expected 'A', 'S', 'V', or 'D', got '{inputs['ptype']}'"
        )

    atomtype = AtomType(**inputs)

    _ensure_unique(
        top.atomtypes,
        lambda at: at.name == atomtype.name and at.ifdef_state == atomtype.ifdef_state,
        f"Duplicate atomtype '{atomtype.name}' found in topology {top.source}.",
    )

    return atomtype


def parse_moleculetype(
    parts: List[str],
    top: Topology,
    ifdef_state: Optional[str] = None
) -> MoleculeType:

    inputs = {
        "name": parts[0],
        "nrexcl": int(parts[1]),
        'ifdef_state': ifdef_state
    }

    moleculetype = MoleculeType(**inputs)

    _ensure_unique(
        top.moleculetypes,
        lambda mt: mt.name == moleculetype.name
        and mt.ifdef_state == moleculetype.ifdef_state,
        f"Duplicate moleculetype '{moleculetype.name}' found in topology {top.source}.",
    )

    return moleculetype


def parse_atom(
    parts: List[str],
    top: Topology,
    mol: MoleculeType,
    ifdef_state: Optional[str] = None
) -> Atom:

    inputs = _parse_atom_values(parts, top)
    inputs["ifdef_state"] = ifdef_state

    if mol.get_atom_by_idx(inputs['nr']) is not None and ifdef_state == "free":
        raise ValueError(
            f"Duplicate atom index '{inputs['nr']}' "
            f"in molecule '{mol.name}' in topology {top.source}"
        )

    atomtype = _find_atomtype(top, inputs['type'])

    inputs["residue"] = mol
    inputs["type"] = atomtype

    return Atom(**inputs)


def parse_exclusion(
    parts: List[str],
    mol: MoleculeType,
    ifdef_state: Optional[str] = None
) -> Exclusion:

    # Resolve atoms
    excluded = _resolve_atoms(parts, mol, "exclusions")
    return Exclusion(excluded, ifdef_state=ifdef_state)


def parse_system(
    parts: List[str],
    top: Topology,
) -> System:

    description = " ".join(parts) if parts else "system"

    if top.system is not None:
        raise ValueError("Multiple [ system ] sections found.")

    system = System(description)
    return system


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
    ifdef_state: Optional[str] = None
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
            f"with atomtypes ({', '.join(a.name for a in atomtypes)})."
        )

    return section_cls(*atomtypes, func, params, ifdef_state=ifdef_state)


def parse_interaction(
    parts: List[str],
    interaction_spec: InteractionSpec,
    section_cls: object,
    top: Topology,
    mol: MoleculeType,
    section: str = "",
    ifdef_state: Optional[str] = None
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
        return [section_cls(*atoms, func, pt.params) for pt in paramtypes]

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

    if section in MOLECULE_SECTION_SPECS and active_molecule is None:
        raise ValueError(
            f"[ {section} ] section found before molecule type in {top.source}."
        )

    if section in PARAMETER_SECTION_SPECS and active_molecule is not None:
        raise ValueError(
            f"[ {section} ] section found inside molecule type "
            f"{active_molecule.name} in {top.source}. "
            f"Parameter types must be defined globally."
        )

    if section == "defaults":
        top.defaults = parse_defaults(parts, top, ifdef_state)
        return active_molecule

    if section == "atomtypes":
        top.atomtypes.append(parse_atomtype(parts, top, ifdef_state))
        return active_molecule

    if section in PARAMETER_SECTION_SPECS:
        interaction_spec, section_cls = PARAMETER_SECTION_SPECS[section]
        interaction_type = parse_interaction_type(
            parts,
            interaction_spec,
            section_cls,
            top,
            section,
            ifdef_state,
        )
        getattr(top, section).append(interaction_type)
        return active_molecule

    if section == "moleculetype":
        molecule = parse_moleculetype(parts, top, ifdef_state)
        top.moleculetypes.append(molecule)
        return molecule

    if section == "atoms":
        active_molecule.atoms.append(
            parse_atom(parts, top, active_molecule, ifdef_state)
        )
        return active_molecule

    if section in MOLECULE_SECTION_SPECS:
        interaction_spec, section_cls = MOLECULE_SECTION_SPECS[section]
        interactions = parse_interaction(
            parts,
            interaction_spec,
            section_cls,
            top,
            active_molecule,
            section,
            ifdef_state,
        )
        getattr(active_molecule, section).extend(interactions)
        return active_molecule

    if section == "exclusions":
        active_molecule.exclusions.append(
            parse_exclusion(parts, active_molecule, ifdef_state)
        )
        return active_molecule

    if section == "system":
        top.system = parse_system(parts, top)
        return active_molecule

    if section == "molecules":
        molecule, count = parse_molecule(parts, top)
        top.molecules[molecule.name] = (molecule, count)
        return active_molecule

    return active_molecule


PARAMETER_SECTION_SPECS: Dict[str, SectionSpec] = {
    "bondtypes": (BONDS, BondType),
    "pairtypes": (PAIRS, PairType),
    "angletypes": (ANGLES, AngleType),
    "dihedraltypes": (DIHEDRALS, DihedralType),
    "constrainttypes": (CONSTRAINTS, ConstraintType),
    "nonbond_params": (NONBOND_PARAMS, NonBondParam),
}

MOLECULE_SECTION_SPECS: Dict[str, SectionSpec] = {
    "bonds": (BONDS, Bond),
    "pairs": (PAIRS, Pair),
    "pairs_nb": (PAIRS_NB, PairNB),
    "angles": (ANGLES, Angle),
    "dihedrals": (DIHEDRALS, Dihedral),
    "constraints": (CONSTRAINTS, Constraint),
    "settles": (SETTLES, Settle),
    "virtual_sites1": (VIRTUAL_SITES1, VirtualSite1),
    "virtual_sites2": (VIRTUAL_SITES2, VirtualSite2),
    "virtual_sites3": (VIRTUAL_SITES3, VirtualSite3),
    "virtual_sites4": (VIRTUAL_SITES4, VirtualSite4),
    "virtual_sitesn": (VIRTUAL_SITESN, VirtualSiteN),
    "dummies1": (VIRTUAL_SITES1, VirtualSite1),
    "dummies2": (VIRTUAL_SITES2, VirtualSite2),
    "dummies3": (VIRTUAL_SITES3, VirtualSite3),
    "dummies4": (VIRTUAL_SITES4, VirtualSite4),
    "dummiesn": (VIRTUAL_SITESN, VirtualSiteN),
    "position_restraints": (POSITION_RESTRAINTS, PositionRestraint),
    "distance_restraints": (DISTANCE_RESTRAINTS, DistanceRestraint),
    "dihedral_restraints": (DIHEDRAL_RESTRAINTS, DihedralRestraint),
    "orientation_restraints": (ORIENTATION_RESTRAINTS, OrientationRestraint),
    "angle_restraints": (ANGLE_RESTRAINTS, AngleRestraint),
    "angle_restraints_z": (ANGLE_RESTRAINTS_Z, AngleRestraintZ),
}
