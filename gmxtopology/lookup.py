from dataclasses import replace
from typing import List, Literal, Mapping, Tuple

from .topology import (
    AngleType,
    Atom,
    BondType,
    ConstraintType,
    DihedralType,
    MoleculeType,
    PairType,
    Topology,
)


def _match_wildcard(pattern: List[str], target: List[str]) -> bool:
    """Check if target matches pattern with 'X' wildcards."""
    return all(p == "X" or p == t for p, t in zip(pattern, target))


def is_exact(dt_names: List[str], at_forward: List[str], at_reverse: List[str]) -> bool:
    """Check for exact match of atom type names."""
    return dt_names == at_forward or dt_names == at_reverse


def is_wild(dt_names: List[str], at_forward: List[str], at_reverse: List[str]) -> bool:
    """Check for wildcard match of atom type names."""
    return _match_wildcard(dt_names, at_forward) or _match_wildcard(
        dt_names,
        at_reverse,
    )


def specificity(dt_names: List[str]) -> int:
    """Calculate specificity of a dihedral type based on non-wildcard entries."""
    return sum(name != "X" for name in dt_names)


Collection = (
    List[BondType]
    | List[PairType]
    | List[AngleType]
    | List[DihedralType]
    | List[ConstraintType]
)


def lookup_paramtype(
    top: Topology,
    *atoms: Atom,
    func: int,
    section: Literal[
        "bonds",
        "pairs",
        "pairs_nb",
        "angles",
        "dihedrals",
        "constraints",
    ],
) -> List[Collection]:
    """Lookup the parameter type(s) for the given atoms in the topology."""

    at_forward = [a.type.name if isinstance(a, Atom) else a for a in atoms]
    at_reverse = list(reversed(at_forward))
    matched_types = None
    if len(atoms) == 2:
        if section in {"pairs", "pairs_nb"}:
            name = "pair"
            for pairtype in top.pairtypes:
                type_names = [pairtype.ai.name, pairtype.aj.name]
                if type_names == at_forward or type_names == at_reverse:
                    matched_types = pairtype
                    break
        elif section == "constraints":
            name = "constraint"
            for constrainttype in top.constrainttypes:
                type_names = [constrainttype.ai.name, constrainttype.aj.name]
                if constrainttype.func != func:
                    continue
                if type_names == at_forward or type_names == at_reverse:
                    matched_types = constrainttype
                    break
        else:
            name = "bond"
            for bondtype in top.bondtypes:
                type_names = [bondtype.ai.name, bondtype.aj.name]
                if bondtype.func != func:
                    continue
                if type_names == at_forward or type_names == at_reverse:
                    matched_types = bondtype
                    break
    elif len(atoms) == 3:
        name = "angle"
        for angletype in top.angletypes:
            type_names = [angletype.ai.name, angletype.aj.name, angletype.ak.name]
            if angletype.func != func:
                continue
            if type_names == at_forward or type_names == at_reverse:
                matched_types = angletype
                break
    elif len(atoms) == 4:
        name = "dihedral"
        candidates = [item for item in top.dihedraltypes if item.func == func]

        exact_matches = []
        wildcard_matches = []
        best_specificity = -1

        for dihedraltype in candidates:
            type_names = [
                dihedraltype.ai.name,
                dihedraltype.aj.name,
                dihedraltype.ak.name,
                dihedraltype.al.name,
            ]

            if is_exact(type_names, at_forward, at_reverse):
                exact_matches.append(dihedraltype)
                continue

            if is_wild(type_names, at_forward, at_reverse):
                matched_specificity = specificity(type_names)
                if matched_specificity > best_specificity:
                    wildcard_matches = [dihedraltype]
                    best_specificity = matched_specificity
                elif matched_specificity == best_specificity:
                    wildcard_matches.append(dihedraltype)

        picked = exact_matches if exact_matches else wildcard_matches

        seen_mult: set[int] = set()
        best: list[DihedralType] = []
        for dihedraltype in picked:
            mult = dihedraltype.params.get("mult")
            if mult is None:
                best.append(dihedraltype)
                continue
            if mult in seen_mult:
                continue
            seen_mult.add(mult)
            best.append(dihedraltype)

        matched_types = best
    else:
        raise ValueError(f"Unsupported interaction size {len(atoms)} in {section}.")

    if not matched_types:
        raise ValueError(
            f"No {name} found for atoms {' '.join(at_forward)} "
            f"with function {func}."
        )

    if not isinstance(matched_types, list):
        matched_types = [matched_types]

    return matched_types


def reduce_atoms(mol: MoleculeType) -> Tuple[List[Atom], Mapping[int, Atom]]:
    """Remove virtual site atoms from the molecule and update all references."""
    new_atoms = []
    old2new = {}
    indices = {
        vsite.ai.nr
        for section in mol._vsite_sections.values()
        for vsite in section
    }

    for atom in mol.atoms:
        if atom.nr in indices:
            continue
        new_nr = len(new_atoms) + 1
        new_atom = replace(atom, nr=new_nr)
        new_atoms.append(new_atom)
        old2new[atom.nr] = new_atom

    return new_atoms, old2new


def reduce_bonded(
    mol: MoleculeType,
    old2new: Mapping[int, Atom],
) -> Mapping[str, List]:
    """Update all connection references in the molecule based on old2new mapping."""
    new_sections = {}
    for section in mol._connection_sections:
        new_list = []
        for item in getattr(mol, section):
            old_ids = [
                getattr(item, attr).nr
                for attr in ("ai", "aj", "ak", "al")
                if hasattr(item, attr)
            ]
            kept_atoms = [old2new[index] for index in old_ids if index in old2new]
            if len(kept_atoms) == len(old_ids):
                for attr, atom in zip(("ai", "aj", "ak", "al"), kept_atoms):
                    if hasattr(item, attr):
                        item = replace(item, **{attr: atom})
                new_list.append(item)
        new_sections[section] = new_list

    return new_sections


def reduce_exclusions(mol: MoleculeType, old2new: Mapping[int, Atom]) -> List:
    """Update exclusions in the molecule based on old2new mapping."""
    new_exclusions = []
    for exclusion in mol.exclusions:
        old_ids = [atom.nr for atom in exclusion.excluded]
        kept_atoms = [old2new[index] for index in old_ids if index in old2new]
        if len(kept_atoms) >= 2:
            new_exclusions.append(replace(exclusion, excluded=tuple(kept_atoms)))
    return new_exclusions


def remove_vsites(mol: MoleculeType) -> None:
    """Remove virtual site atoms from the molecule and update all references."""
    new_atoms, old2new = reduce_atoms(mol)
    new_bonded = reduce_bonded(mol, old2new)
    new_exclusions = reduce_exclusions(mol, old2new)

    mol.atoms = new_atoms
    for field_name, new_list in new_bonded.items():
        setattr(mol, field_name, new_list)
    mol.exclusions = new_exclusions

    for vsite_kind in mol._vsite_sections:
        setattr(mol, vsite_kind, [])
