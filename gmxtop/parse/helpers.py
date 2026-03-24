from typing import List, Tuple, Mapping, Literal
from dataclasses import replace
from ..topology import (
    Atom, BondType, PairType, AngleType, DihedralType, ConstraintType,
    MoleculeType, Topology
)


def _match_wildcard(pattern: List[str], target: List[str]) -> bool:
    """Check if target matches pattern with 'X' wildcards."""
    return all(p == "X" or p == t for p, t in zip(pattern, target))


def is_exact(dt_names: List[str], at_forward: List[str], at_reverse: List[str]) -> bool:
    """Check for exact match of atom type names."""
    return dt_names == at_forward or dt_names == at_reverse


def is_wild(dt_names: List[str], at_forward: List[str], at_reverse: List[str]) -> bool:
    """Check for wildcard match of atom type names."""
    return (
        _match_wildcard(dt_names, at_forward) or _match_wildcard(dt_names, at_reverse)
    )


def specificity(dt_names: List[str]) -> int:
    """Calculate specificity of a dihedral type based on non-wildcard entries."""
    return sum(x != "X" for x in dt_names)


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
            for pt in top.pairtypes:
                pt_list = [pt.ai.name, pt.aj.name]
                if pt_list == at_forward or pt_list == at_reverse:
                    matched_types = pt
                    break
        elif section == "constraints":
            name = "constraint"
            for ct in top.constrainttypes:
                ct_list = [ct.ai.name, ct.aj.name]
                if ct.func != func:
                    continue
                if ct_list == at_forward or ct_list == at_reverse:
                    matched_types = ct
                    break
        else:
            name = "bond"
            for bt in top.bondtypes:
                bt_list = [bt.ai.name, bt.aj.name]
                if bt.func != func:
                    continue
                if bt_list == at_forward or bt_list == at_reverse:
                    matched_types = bt
                    break

    elif len(atoms) == 3:
        name = "angle"
        for at in top.angletypes:
            at_list = [at.ai.name, at.aj.name, at.ak.name]
            if at.func != func:
                continue
            if at_list == at_forward or at_list == at_reverse:
                matched_types = at
                break

    elif len(atoms) == 4:
        name = "dihedral"

        # Only consider dihedrals with matching func
        candidates = [dt for dt in top.dihedraltypes if dt.func == func]

        # Prefer exact matches; otherwise use wildcard matches with max specificity
        exact = []
        wild = []
        best_spec = -1

        for dt in candidates:
            dt_names = [dt.ai.name, dt.aj.name, dt.ak.name, dt.al.name]

            if is_exact(dt_names, at_forward, at_reverse):
                exact.append(dt)
                continue

            if is_wild(dt_names, at_forward, at_reverse):
                spec = specificity(dt_names)
                if spec > best_spec:
                    wild = [dt]
                    best_spec = spec
                elif spec == best_spec:
                    wild.append(dt)

        picked = exact if exact else wild

        # Dedupe: keep at most one term per multiplicity (mult); allow mult=None entries
        seen_mult: set[int] = set()
        best: list[DihedralType] = []
        for dt in picked:
            mult = dt.params.get("mult")
            if mult is None:
                best.append(dt)
                continue
            if mult in seen_mult:
                continue
            seen_mult.add(mult)
            best.append(dt)

        matched_types = best

    if not matched_types:
        raise ValueError(
            f"No {name} found for atoms "
            f"{' '.join(at_forward)} with function {func}."
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
    old2new: Mapping[int, Atom]
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
            kept_atoms = [old2new[i] for i in old_ids if i in old2new]
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
    for ex in mol.exclusions:
        old_ids = [a.nr for a in ex.excluded]
        kept_atoms = [old2new[i] for i in old_ids if i in old2new]
        if len(kept_atoms) >= 2:
            new_exclusions.append(replace(ex, excluded=tuple(kept_atoms)))
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
