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

This installs `ipykernel` as well, so the environment is ready for Jupyter and
IPython kernel use.

From source:

```bash
pip install git+https://github.com/vojtechkostal/gmxtopology.git
```

## Usage
See [examples/example.ipynb](examples/example.ipynb) for a fuller walkthrough.

```python
from gmxtopology import Topology

# load topology from file
fn_top = "./examples/topol.top"
top = Topology(fn_top)

# modifying topology parameters
for atom in top.molecules['POPC'][0].atoms:
    atom.update(charge=0.0)

for atomtype in top.atomtypes:
    atomtype.update(sigma=0.31)

# remove virtual sites from a molecule
top.molecules['SOL'][0].remove_vsites()

# write topology into a file
top.write("./examples/topol-processed.top", overwrite=True)
```

## Publishing
GitHub Trusted Publishing is set up via
[`publish.yml`](.github/workflows/publish.yml). The one-time PyPI and TestPyPI
configuration steps are documented in [docs/publishing.md](docs/publishing.md).
