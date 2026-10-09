# Changelog

## 0.3.0

Breaking release: the reader, data model and writer were rewritten as a linear
pipeline (preprocess, parse, resolve, write).

### Migration

- `top.molecules` is now an ordered list of `(MoleculeType, count)`; repeated
  molecule names are kept. Look up a molecule type with
  `top.moleculetype("NAME")`.
- Conditional blocks are evaluated: pass `defines={"NAME": None}` to
  `Topology`. The `ifdef_state` of records and the round-tripping of
  `#ifdef` blocks are gone.
- The per-directive record classes (`Bond`, `Angle`, ...) are replaced by
  `Interaction`; parameter tables hold `TypeEntry` records. Atoms and
  interactions refer to atom types by name and to atoms by number, e.g.
  `atom.type` is a string and `bond.atoms == (1, 2)`.
- `molecule.remove_vsites()` is now `drop_vsites(molecule)`.
- `RawSection` and `raw_sections` are gone: CMAP is parsed and unknown
  directives are errors.

### Fixes

- Bonded parameters are matched on bonded atom types (OPLS-AA, GROMOS).
- Explicit `[ pairtypes ]` are used even with `gen-pairs = yes`.
- Repeated molecule names in `[ molecules ]` are kept in order.
- Residue names are written as read.
- Only the CMAP grids that are used are written, and `[ cmap ]` entries are
  renumbered when virtual sites are removed.
- `[ dihedraltypes ]`: two-type form, GROMACS wildcard matching, and adjacent
  type-9 terms with equal multiplicity.
- Parameters are written with full precision.
- Includes are found through `$GMXLIB` and the GROMACS data directory; errors
  report file and line; `#undef`, line continuation, `[ intermolecular_interactions ]`
  and bonds without parameters are supported.
- Fourier dihedrals (type 5) take four parameters.

## 0.2.0

### Highlights

- Preserve `#define`, `#ifdef`, `#ifndef`, `#else`, and nested conditional
  blocks while reading and writing topology files.
- Support marker defines such as `#define POSRES` and replacement text with
  multiple tokens.
- Parse optional and variable-width GROMACS records, including topology-B
  free-energy parameters and additional virtual-site forms.
- Preserve unsupported force-field-specific sections as opaque records,
  including CHARMM CMAP tables.
- Write flattened molecular interactions with resolved bonded parameters,
  filtered atom types, and relevant nonbond overrides.
- Match GROMACS molecule names case-insensitively.
- Speed up large CHARMM-derived topology parsing with lazy atom-type and
  bonded-parameter indexes. The realistic prosECCo75 template used during
  development dropped from roughly 33 seconds to under 1 second locally.
- Add regression fixtures from official GROMACS files, prosECCo75, and the
  ParmEd `GmxTests` corpus. Flattened fixtures are validated with native
  GROMACS during release preparation.
- Keep PyPI artifacts lightweight: distributions include one self-contained
  notebook topology while larger examples and regression fixtures remain
  available from the repository.

## 0.1.1

- Initial PyPI release under the `gmxtopology` package name.
