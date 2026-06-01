# gmxtopology

## Description
`gmxtopology` is a Python package for parsing and editing GROMACS topology
files.

## Package layout
The main flow is organized around a few top-level modules:

- `gmxtopology.topology`: topology data models
- `gmxtopology.io`: topology file reading and writing
- `gmxtopology.parser`: section parsing and dispatch
- `gmxtopology.lookup`: parameter lookup and virtual-site reduction helpers
- `gmxtopology.schema`: schema objects used by the parser
- `gmxtopology.interaction_specs`: GROMACS interaction tables

## Installation
```bash
pip install gmxtopology
```

For Jupyter/IPython kernel support, install the optional notebook extra:

```bash
pip install "gmxtopology[notebook]"
```

From source:

```bash
pip install git+https://github.com/vojtechkostal/gmxtopology.git
```

## Usage
See [examples/example.ipynb](examples/example.ipynb) for a fuller walkthrough.
The packaged walkthrough uses the single self-contained
[examples/example.top](examples/example.top) file.

```python
from gmxtopology import Topology

# load topology from file
fn_top = "./examples/example.top"
top = Topology(fn_top)

# modifying topology parameters
for atom in top.molecules['MOL'][0].atoms:
    atom.update(charge=0.0)

for atomtype in top.atomtypes:
    atomtype.update(sigma=0.31)

# remove virtual sites from a molecule
top.molecules['SOL'][0].remove_vsites()

# write topology into a file
top.write("./examples/topol-processed.top", overwrite=True)
```

## Topology preprocessing
The reader follows local `#include` directives relative to the including file
and preserves conditional topology sections when writing a modified topology.
Supported directives include:

- `#define NAME` and `#define NAME replacement text`
- `#ifdef NAME` and `#ifndef NAME`
- `#else` and `#endif`, including nested blocks
- `#include "relative/path.itp"`

Missing includes remain errors unless they occur inside a conditional block.
This allows topologies with absent optional restraint files to load when those
files are referenced only behind flags such as `#ifdef REST_ON`.

For free-energy topologies, `[ atoms ]` records and supported bonded
interactions may also include the optional topology-B parameters documented by
GROMACS. They are exposed with `_b` suffixes, such as `charge_b`, `b0_b`, and
`kb_b`.

Sections that gmxtopology does not edit are carried through as opaque records
when writing a flattened topology. This preserves force-field-specific data
such as CHARMM `[ cmaptypes ]`, `[ implicit_genborn_params ]`, and molecular
`[ cmap ]` sections.

Written topologies contain only the `[ atomtypes ]` used by molecules in the
`[ molecules ]` section or referenced by preserved opaque global sections.
Bonded parameter tables such as `[ bondtypes ]`, `[ angletypes ]`, and
`[ dihedraltypes ]` are resolved onto molecular interactions and omitted from
the output. Relevant `[ nonbond_params ]` entries remain global because they
cannot be expanded onto bonded records.

## Testing
Run the parser regression suite with:

```bash
python -m pytest
```

The repository test suite includes pinned
[official GROMACS fixtures](tests/fixtures/gromacs-v2026.2/SOURCE.md) covering
SPC/E, TIP3P, TIP4P, and urea topologies, plus
[prosECCo75 CHARMM36-derived fixtures](tests/fixtures/prosECCo75-e4831a4/SOURCE.md)
covering POPC and CMAP tables. A curated
[ParmEd `GmxTests` snapshot](tests/fixtures/parmed-96ec61a/SOURCE.md) adds
real-world water, peptide, solvated protein, mixed-solvent, and DPPC cases.

PyPI distributions intentionally contain only the package code, documentation,
the example notebook, and its small self-contained topology. Clone the
repository explicitly when you need the larger examples and regression
fixtures:

```bash
git clone --depth 1 https://github.com/vojtechkostal/gmxtopology.git
```

## Changelog
See [CHANGELOG.md](CHANGELOG.md) for release highlights.

## Publishing
GitHub Trusted Publishing is set up via
[`publish.yml`](.github/workflows/publish.yml). The one-time PyPI and TestPyPI
configuration steps are documented in [docs/publishing.md](docs/publishing.md).
