# gmxtop

## Description
`gmxtop` is a Python package for parsing and editing GROMACS topology files.

## Package layout
The main flow is organized around a few top-level modules:

- `gmxtop.topology`: topology data models
- `gmxtop.io`: topology file reading and writing
- `gmxtop.parser`: parser helpers and section parsers
- `gmxtop.schema`: schema objects used by the parser
- `gmxtop.interaction_specs`: GROMACS interaction tables

## Installation
```bash
pip install gmxtop
```

From source:

```bash
pip install git+https://github.com/vojtechkostal/gmxtop.git
```

## Usage
See [examples/example.ipynb](/home/vojta/code/gmxtop-git/examples/example.ipynb) for a fuller walkthrough.

```python
from gmxtop import Topology

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
