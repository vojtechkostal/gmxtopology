# ParmEd GmxTests topology fixtures

These files are a deliberately small snapshot of the
[ParmEd `GmxTests` corpus](https://github.com/ParmEd/ParmEd/tree/master/GmxTests)
at commit
[`96ec61a3c83be5eab2058e3f4e6e24053c8a219d`](https://github.com/ParmEd/ParmEd/tree/96ec61a3c83be5eab2058e3f4e6e24053c8a219d/GmxTests).
They provide real-world GROMACS topologies for parser regression tests.

The copied cases are:

- `01.1water`: TIP3P water and an `#ifndef FLEXIBLE` block.
- `02.6water`: TIP4P water with a virtual site.
- `03.AlaGlu`: a two-peptide system with included parameter and molecule files.
- `07.DHFR-Liquid-NoPBC`: a larger solvated DHFR topology with a local Amber03
  force field.
- `11a.Toluene-Cyclohexane_conversion`: a compact ParmEd-generated mixed
  solvent topology.
- `12.DPPC`: a DPPC bilayer topology with nonbond overrides.
- `12A.DPPC_Amber`: an Amber-style DPPC bilayer topology with a conditional
  improper torsion section.

The substantially larger fully expanded solvent topologies are intentionally
not copied into the test suite. They are useful for occasional stress tests,
but duplicate coverage while adding several megabytes and around 20 seconds
per parse.

ParmEd is distributed under the GNU Lesser General Public License. The
upstream license text is included in [`GNU_LGPL_v2`](GNU_LGPL_v2).
