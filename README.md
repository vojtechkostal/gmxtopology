# gmxtopology

## Description
`gmxtopology` is a Python package for parsing and editing GROMACS topology
files.

## Package layout
A topology is read in four small steps, each in its own module:

- `gmxtopology.preprocess`: the GROMACS preprocessor (`#include`, `#define`,
  `#ifdef`, comments, continuation lines), yielding the active lines
- `gmxtopology.parser`: directive parsing; `gmxtopology.specs` lists the
  parameters of every interaction function type
- `gmxtopology.lookup`: default parameters from `[ bondtypes ]` and friends
- `gmxtopology.model`: the topology data model
- `gmxtopology.reduce` and `gmxtopology.writer`: select what the system needs and
  write it

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
from gmxtopology import Topology, drop_vsites

# load topology from file; defines act like `define = -DPOSRES` in an mdp file
top = Topology("./examples/example.top", defines={"POSRES": None})

# modifying topology parameters
for atom in top.moleculetype("MOL").atoms:
    atom.update(charge=0.0)

for atomtype in top.atomtypes:
    atomtype.update(sigma=0.31)

# remove virtual sites from a molecule
drop_vsites(top.moleculetype("SOL"))

# write topology into a file
top.write("./examples/topol-processed.top", overwrite=True)
```

## Reading
The reader behaves like `grompp`:

- `#include` files are searched next to the including file, in `include_dirs`,
  and in `$GMXLIB` and the GROMACS `share/gromacs/top` directory (found through
  `$GMXDATA` or the `gmx` executable), so `#include "oplsaa.ff/forcefield.itp"`
  works.
- `#ifdef`, `#ifndef`, `#else` and `#endif` are evaluated against the macros
  given as `defines` and the `#define` lines of the files, so only the active
  branches are read. `#define` macros that stand for one number stay symbolic
  in the parameters and can be edited through `top.defines`; macros that stand
  for several numbers (such as GROMOS `gb_21`) are expanded.
- Parameters left out of `[ bonds ]`, `[ pairs ]`, `[ angles ]`, `[ dihedrals ]`
  and `[ constraints ]` are filled in from the parameter tables, matching on
  *bonded* atom types, either orientation, `X` wildcards in `[ dihedraltypes ]`
  (including the two-type form), multi-term type-9 dihedrals, and
  `gen-pairs`.
- `[ cmaptypes ]` and `[ cmap ]`, `[ intermolecular_interactions ]`, B-state
  parameters and the other interaction directives of the GROMACS manual are
  parsed into typed records. Unknown directives are errors with file and line.

Molecular interactions are `Interaction(section, atoms, func, params)` records
in `molecule.interactions[section]`, also reachable as `molecule.bonds`,
`molecule.dihedrals`, and so on. Atoms are referred to by their 1-based number
and atom types by name. `top.molecules` is the ordered `[ molecules ]` list of
`(MoleculeType, count)` pairs.

## Writing
Written topologies are flat and contain only what the simulated system needs:

- the molecule types listed in `[ molecules ]` (in that order),
- the atom types used by their atoms, and the `[ nonbond_params ]` between them,
- the CMAP grids that a `[ cmap ]` of those molecules refers to,
- macros that parameters still refer to.

Bonded parameter tables are not written, because the parameters were resolved
onto the molecular interactions. Includes and conditional blocks do not appear.

## Testing
Run the parser regression suite with:

```bash
python -m pytest
```

Tests that compare a topology with its written version through `gmx grompp`
and `gmx dump` run when `gmx` is on the `PATH`.

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
