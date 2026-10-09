"""Which parts of a topology matter for the simulated system."""

from dataclasses import replace

from .model import AtomType, CmapType, MoleculeType, Topology, TypeEntry
from .specs import VSITE_SECTIONS


def used_moleculetypes(top: Topology) -> list[MoleculeType]:
    """The molecule types listed in ``[ molecules ]``, in order of first use."""
    used: list[MoleculeType] = []
    for mol, _ in top.molecules:
        if not any(mol is other for other in used):
            used.append(mol)
    return used


def used_atomtypes(top: Topology) -> list[AtomType]:
    """Atom types of the atoms of the system (A and B state)."""
    names = {
        name
        for mol in used_moleculetypes(top)
        for atom in mol.atoms
        for name in (atom.type, atom.type_b)
        if name is not None
    }
    return [atomtype for atomtype in top.atomtypes if atomtype.name in names]


def used_nonbond_params(top: Topology, atomtypes: list[AtomType]) -> list[TypeEntry]:
    """``[ nonbond_params ]`` overrides between atom types that are used."""
    names = {atomtype.name for atomtype in atomtypes}
    return [
        entry
        for entry in top.types.get("nonbond_params", [])
        if all(name in names for name in entry.types)
    ]


def used_cmaptypes(top: Topology) -> list[CmapType]:
    """CMAP grids that the ``[ cmap ]`` entries of the system refer to."""
    btypes = {atomtype.name: atomtype.btype for atomtype in top.atomtypes}
    wanted: set[tuple[str, ...]] = set()
    for mol in used_moleculetypes(top):
        for cmap in mol.interactions.get("cmap", []):
            types = tuple(btypes[mol.atoms[nr - 1].type] for nr in cmap.atoms)
            wanted.update((types, types[::-1]))
    return [cmaptype for cmaptype in top.cmaptypes if cmaptype.types in wanted]


def virtual_site_numbers(mol: MoleculeType) -> set[int]:
    """Numbers of the atoms of ``mol`` that are virtual sites."""
    return {
        interaction.atoms[0]
        for section in VSITE_SECTIONS
        for interaction in mol.interactions.get(section, [])
    }


def drop_vsites(mol: MoleculeType) -> None:
    """Remove virtual sites from a molecule and renumber what remains.

    Interactions that involve a removed atom are dropped.
    """
    sites = virtual_site_numbers(mol)
    new_number: dict[int, int] = {}
    atoms = []
    for atom in mol.atoms:
        if atom.nr not in sites:
            new_number[atom.nr] = len(atoms) + 1
            atoms.append(replace(atom, nr=len(atoms) + 1))

    interactions = {}
    for section, items in mol.interactions.items():
        if section in VSITE_SECTIONS:
            continue
        interactions[section] = [
            replace(item, atoms=tuple(new_number[nr] for nr in item.atoms))
            for item in items
            if all(nr in new_number for nr in item.atoms)
        ]
    exclusions = [
        kept
        for kept in (
            tuple(new_number[nr] for nr in group if nr in new_number)
            for group in mol.exclusions
        )
        if len(kept) >= 2
    ]
    mol.atoms, mol.interactions, mol.exclusions = atoms, interactions, exclusions
