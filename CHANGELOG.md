# Changelog

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
