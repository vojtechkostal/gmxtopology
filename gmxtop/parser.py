"""Single entrypoint for topology parsing and parser-side transformations.

This module gathers the parsing functions and parser-side transformations so
the main read flow is easier to follow from one place.
"""

from .parse.helpers import (
    is_exact,
    is_wild,
    lookup_paramtype,
    reduce_atoms,
    reduce_bonded,
    reduce_exclusions,
    remove_vsites,
    specificity,
)
from .parse.parsers import (
    MOLECULE_SECTION_SPECS,
    PARAMETER_SECTION_SPECS,
    apply_section_line,
    parse_atom,
    parse_atomtype,
    parse_defaults,
    parse_exclusion,
    parse_interaction,
    parse_interaction_type,
    parse_molecule,
    parse_moleculetype,
    parse_system,
)

__all__ = [
    "MOLECULE_SECTION_SPECS",
    "PARAMETER_SECTION_SPECS",
    "apply_section_line",
    "is_exact",
    "is_wild",
    "lookup_paramtype",
    "parse_atom",
    "parse_atomtype",
    "parse_defaults",
    "parse_exclusion",
    "parse_interaction",
    "parse_interaction_type",
    "parse_molecule",
    "parse_moleculetype",
    "parse_system",
    "reduce_atoms",
    "reduce_bonded",
    "reduce_exclusions",
    "remove_vsites",
    "specificity",
]
