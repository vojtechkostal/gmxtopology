# Official GROMACS topology fixtures

These files are a deliberately small snapshot of the official GROMACS topology
distribution. They were downloaded from the
[GROMACS `v2026.2` tag](https://gitlab.com/gromacs/gromacs/-/tree/v2026.2/share/top/amber99.ff)
and are used to test gmxtopology against real topology files accepted by
GROMACS.

The copied upstream files are:

- [`amber99.ff/forcefield.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/forcefield.itp)
- [`amber99.ff/ffnonbonded.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/ffnonbonded.itp)
- [`amber99.ff/ffbonded.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/ffbonded.itp)
- [`amber99.ff/spce.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/spce.itp)
- [`amber99.ff/tip3p.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/tip3p.itp)
- [`amber99.ff/tip4p.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/tip4p.itp)
- [`amber99.ff/urea.itp`](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/share/top/amber99.ff/urea.itp)

The `topol-*.top`, `*.gro`, and `grompp.mdp` files in this directory are
minimal gmxtopology test wrappers. The upstream GROMACS license is available
in the local [`COPYING`](COPYING) file and the
[official source tree](https://gitlab.com/gromacs/gromacs/-/blob/v2026.2/COPYING).
