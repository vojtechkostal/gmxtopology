from dataclasses import dataclass, replace
from typing import Literal, Mapping

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

ATOM_FIELDS = ("ai", "aj", "ak", "al")
PARAMTYPE_SECTIONS = {
    "bonds": ("bond", "bondtypes"),
    "pairs": ("pair", "pairtypes"),
    "pairs_nb": ("pair", "pairtypes"),
    "angles": ("angle", "angletypes"),
    "constraints": ("constraint", "constrainttypes"),
}
ParamType = BondType | PairType | AngleType | DihedralType | ConstraintType


@dataclass
class ParamTypeIndex:
    size: int
    exact: dict[tuple[int, tuple[str, ...]], list[ParamType]]
    wildcards: dict[int, list[DihedralType]]


def _atomtype_names(paramtype: ParamType, n_atoms: int) -> tuple[str, ...]:
    return tuple(
        getattr(paramtype, field).name
        for field in ATOM_FIELDS[:n_atoms]
    )


def _canonical_atomtypes(names: tuple[str, ...]) -> tuple[str, ...]:
    reverse = tuple(reversed(names))
    return min(names, reverse)


def _matches_wildcard(pattern: tuple[str, ...], target: tuple[str, ...]) -> bool:
    return all(
        expected == "X" or expected == actual
        for expected, actual in zip(pattern, target)
    )


def _paramtype_index(
    top: Topology,
    collection_name: str,
    n_atoms: int,
) -> ParamTypeIndex:
    collection = getattr(top, collection_name)
    indexes = getattr(top, "_paramtype_indexes", {})
    cached = indexes.get(collection_name)
    if cached is not None and cached.size == len(collection):
        return cached

    exact = {}
    wildcards = {}
    for paramtype in collection:
        names = _atomtype_names(paramtype, n_atoms)
        if "X" in names:
            wildcards.setdefault(paramtype.func, []).append(paramtype)
        else:
            key = paramtype.func, _canonical_atomtypes(names)
            exact.setdefault(key, []).append(paramtype)

    cached = ParamTypeIndex(
        size=len(collection),
        exact=exact,
        wildcards=wildcards,
    )
    indexes[collection_name] = cached
    top._paramtype_indexes = indexes
    return cached


def _lookup_dihedraltypes(
    top: Topology,
    at_forward: tuple[str, ...],
    func: int,
) -> list[DihedralType]:
    index = _paramtype_index(top, "dihedraltypes", 4)
    exact_matches = index.exact.get((func, _canonical_atomtypes(at_forward)), [])
    wildcard_matches = []
    best_specificity = -1
    at_reverse = tuple(reversed(at_forward))

    for dihedraltype in index.wildcards.get(func, []):
        type_names = _atomtype_names(dihedraltype, 4)
        if not _matches_wildcard(type_names, at_forward):
            if not _matches_wildcard(type_names, at_reverse):
                continue

        matched_specificity = sum(name != "X" for name in type_names)
        if matched_specificity > best_specificity:
            wildcard_matches = [dihedraltype]
            best_specificity = matched_specificity
        elif matched_specificity == best_specificity:
            wildcard_matches.append(dihedraltype)

    matched_types = exact_matches or wildcard_matches
    unique_types = []
    seen_multiplicities = set()
    for dihedraltype in matched_types:
        multiplicity = dihedraltype.params.get("mult")
        if multiplicity is not None and multiplicity in seen_multiplicities:
            continue
        if multiplicity is not None:
            seen_multiplicities.add(multiplicity)
        unique_types.append(dihedraltype)
    return unique_types


InteractionSection = Literal[
    "bonds",
    "pairs",
    "pairs_nb",
    "angles",
    "dihedrals",
    "constraints",
]


def lookup_paramtype(
    top: Topology,
    *atoms: Atom,
    func: int,
    section: InteractionSection,
) -> list[ParamType]:
    """Look up global parameters for a molecular interaction."""

    at_forward = tuple(atom.type.name for atom in atoms)

    if section == "dihedrals":
        name = "dihedral"
        matched_types = _lookup_dihedraltypes(top, at_forward, func)
    else:
        try:
            name, collection_name = PARAMTYPE_SECTIONS[section]
        except KeyError as exc:
            raise ValueError(
                f"Unsupported interaction size {len(atoms)} in {section}."
            ) from exc

        index = _paramtype_index(top, collection_name, len(atoms))
        matched_types = index.exact.get(
            (func, _canonical_atomtypes(at_forward)),
            [],
        )[:1]

    if not matched_types:
        raise ValueError(
            f"No {name} found for atoms {' '.join(at_forward)} "
            f"with function {func}."
        )
    return matched_types


def reduce_atoms(mol: MoleculeType) -> tuple[list[Atom], Mapping[int, Atom]]:
    """Remove virtual site atoms from the molecule and update all references."""
    new_atoms = []
    old2new = {}
    indices = {
        vsite.ai.nr
        for section in mol.VSITE_SECTIONS
        for vsite in getattr(mol, section)
    }

    for atom in mol.atoms:
        if atom.nr in indices:
            continue
        new_atom = replace(atom, nr=len(new_atoms) + 1)
        new_atoms.append(new_atom)
        old2new[atom.nr] = new_atom

    return new_atoms, old2new


def reduce_connections(
    mol: MoleculeType,
    old2new: Mapping[int, Atom],
) -> dict[str, list]:
    """Update all connection references in the molecule based on old2new mapping."""
    new_sections = {}
    for section in mol.CONNECTION_SECTIONS:
        new_items = []
        for item in getattr(mol, section):
            atom_fields = [field for field in ATOM_FIELDS if hasattr(item, field)]
            atoms = {
                field: old2new[getattr(item, field).nr]
                for field in atom_fields
                if getattr(item, field).nr in old2new
            }
            if len(atoms) == len(atom_fields):
                new_items.append(replace(item, **atoms))
        new_sections[section] = new_items

    return new_sections


def reduce_exclusions(
    mol: MoleculeType,
    old2new: Mapping[int, Atom],
) -> list:
    """Update exclusions in the molecule based on old2new mapping."""
    new_exclusions = []
    for exclusion in mol.exclusions:
        kept_atoms = [
            old2new[atom.nr]
            for atom in exclusion.excluded
            if atom.nr in old2new
        ]
        if len(kept_atoms) >= 2:
            new_exclusions.append(replace(exclusion, excluded=kept_atoms))
    return new_exclusions


def remove_vsites(mol: MoleculeType) -> None:
    """Remove virtual site atoms from the molecule and update all references."""
    new_atoms, old2new = reduce_atoms(mol)
    new_connections = reduce_connections(mol, old2new)
    new_exclusions = reduce_exclusions(mol, old2new)

    mol.atoms = new_atoms
    for field_name, new_list in new_connections.items():
        setattr(mol, field_name, new_list)
    mol.exclusions = new_exclusions

    for vsite_kind in mol.VSITE_SECTIONS:
        setattr(mol, vsite_kind, [])
