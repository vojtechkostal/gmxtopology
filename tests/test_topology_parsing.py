import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from gmxtopology import Topology
from gmxtopology.topology import Atom, AtomType, Defaults, MoleculeType, System


ROOT = Path(__file__).resolve().parents[1]
GROMACS_FIXTURES = ROOT / "tests/fixtures/gromacs-v2026.2"
PARMED_FIXTURES = ROOT / "tests/fixtures/parmed-96ec61a"
PROSECCO_FIXTURES = ROOT / "tests/fixtures/prosECCo75-e4831a4"


class TopologyParsingTests(unittest.TestCase):
    def test_section_from_line_helpers_construct_records(self) -> None:
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "from-line.top"
            path.write_text("")
            top = Topology.__new__(Topology)
            top.source = path.resolve()
            top.defaults = None
            top.atomtypes = []
            top.nonbond_params = []
            top.bondtypes = []
            top.pairtypes = []
            top.angletypes = []
            top.dihedraltypes = []
            top.constrainttypes = []
            top.moleculetypes = []
            top.system = None
            top.molecules = {}
            top.defines = []

            defaults = Defaults.from_line(["1", "2"], top)
            self.assertEqual(defaults.gen_pairs, "no")

            top.defaults = defaults
            atomtype = AtomType.from_line(
                ["OW", "15.9994", "-0.834", "A", "0.3165", "0.65"],
                top,
            )
            top.atomtypes.append(atomtype)
            self.assertEqual(atomtype.name, "OW")

            molecule = MoleculeType.from_line(["SOL", "3"], top)
            atom = Atom.from_line(
                ["1", "OW", "1", "SOL", "OW", "1", "-0.834", "15.9994"],
                top,
                molecule,
            )
            system = System.from_line(["Water"], top)

        self.assertEqual(molecule.name, "SOL")
        self.assertEqual(atom.type.name, "OW")
        self.assertEqual(system.description, "Water")

    def test_example_topology_parses(self) -> None:
        top = Topology(ROOT / "examples/topol.top")

        self.assertEqual(len(top.atomtypes), 49)
        self.assertEqual(len(top.moleculetypes), 2)
        self.assertEqual(len(top.molecules), 2)
        self.assertEqual(len(top.atoms), 53444)

    def test_example_topology_with_conditionals_parses(self) -> None:
        top = Topology(ROOT / "examples/topol_new.top")

        self.assertEqual(len(top.atomtypes), 22)
        self.assertEqual(len(top.moleculetypes), 2)
        self.assertEqual(len(top.molecules), 2)
        self.assertEqual(len(top.atoms), 53444)

    def test_roundtrip_output_reparses(self) -> None:
        source = Topology(ROOT / "examples/topol.top")
        output = Path("/tmp/gmxtop-roundtrip.top")

        source.write(output, overwrite=True)
        reparsed = Topology(output)

        self.assertEqual(len(reparsed.moleculetypes), len(source.moleculetypes))
        self.assertEqual(len(reparsed.molecules), len(source.molecules))
        self.assertEqual(len(reparsed.atoms), len(source.atoms))
        self.assertEqual(
            {atom.type.name for atom in reparsed.atoms},
            {atom.type.name for atom in source.atoms},
        )

    def test_official_gromacs_water_topologies_roundtrip(self) -> None:
        expected = {
            "spce": (3, 1, 0),
            "tip3p": (3, 1, 0),
            "tip4p": (4, 1, 1),
        }

        with TemporaryDirectory() as tmpdir:
            for model, (atoms, settles, virtual_sites3) in expected.items():
                with self.subTest(model=model):
                    source = GROMACS_FIXTURES / f"topol-{model}.top"
                    output = Path(tmpdir) / f"topol-{model}.top"
                    top = Topology(source)
                    molecule = top.molecules["SOL"][0]

                    self.assertEqual(len(molecule.atoms), atoms)
                    self.assertEqual(len(molecule.settles), settles)
                    self.assertEqual(len(molecule.virtual_sites3), virtual_sites3)

                    top.write(output)
                    reparsed = Topology(output)
                    self.assertEqual(len(reparsed.atoms), atoms)

    def test_official_gromacs_urea_topology_roundtrip(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "topol-urea.top"
            top = Topology(GROMACS_FIXTURES / "topol-urea.top")
            molecule = top.molecules["URE"][0]

            self.assertEqual(len(molecule.atoms), 8)
            self.assertEqual(len(molecule.bonds), 7)
            self.assertTrue(all(bond.func == 1 for bond in molecule.bonds))
            self.assertEqual(len(molecule.dihedrals), 15)

            top.write(output)
            reparsed = Topology(output)
            self.assertEqual(len(reparsed.atoms), 8)

    def test_parmed_real_world_topologies_roundtrip(self) -> None:
        expected = {
            "01.1water/topol.top": 3,
            "02.6water/topol.top": 24,
            "03.AlaGlu/topol.top": 49,
            (
                "11a.Toluene-Cyclohexane_conversion/"
                "toluene_cyclohexane_10_500_parmed.top"
            ): 9150,
            "12.DPPC/topol.top": 1132,
            "12A.DPPC_Amber/topol.top": 1132,
        }

        with TemporaryDirectory() as tmpdir:
            for index, (relative, atoms) in enumerate(expected.items()):
                with self.subTest(topology=relative):
                    top = Topology(PARMED_FIXTURES / relative)
                    output = Path(tmpdir) / f"parmed-{index}.top"

                    self.assertEqual(len(top.atoms), atoms)
                    top.write(output)

                    reparsed = Topology(output)
                    self.assertEqual(len(reparsed.atoms), atoms)

    def test_parmed_solvated_dhfr_topology_parses(self) -> None:
        top = Topology(PARMED_FIXTURES / "07.DHFR-Liquid-NoPBC/topol.top")

        self.assertEqual(len(top.moleculetypes), 12)
        self.assertEqual(len(top.molecules), 3)
        self.assertEqual(len(top.atoms), 23569)

    def test_remove_virtual_sites_updates_atom_references(self) -> None:
        top = Topology(GROMACS_FIXTURES / "topol-tip4p.top")
        molecule = top.molecules["SOL"][0]

        molecule.remove_vsites()

        self.assertEqual(len(molecule.atoms), 3)
        self.assertEqual(len(molecule.virtual_sites3), 0)
        self.assertTrue(
            all(
                atom.nr <= 3
                for exclusion in molecule.exclusions
                for atom in exclusion.excluded
            )
        )

    def test_prosecco_popc_topology_parses(self) -> None:
        top = Topology(PROSECCO_FIXTURES / "topol-popc.top")
        molecule = top.molecules["POPC_s"][0]

        self.assertEqual(len(molecule.atoms), 134)
        self.assertEqual(len(molecule.bonds), 133)
        self.assertEqual(len(molecule.pairs), 356)
        self.assertEqual(len(molecule.angles), 256)
        self.assertEqual(len(molecule.dihedrals), 467)

    def test_prosecco_charmm_sections_roundtrip(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "topol-cmap.top"
            top = Topology(PROSECCO_FIXTURES / "topol-cmap.top")

            self.assertEqual(
                [section.name for section in top.raw_sections],
                ["implicit_genborn_params", "cmaptypes", "cmaptypes"],
            )
            self.assertEqual(
                [section.name for section in top.molecules["CMAP"][0].raw_sections],
                ["cmap"],
            )

            top.write(output)
            reparsed = Topology(output)

        self.assertEqual(
            [section.name for section in reparsed.raw_sections],
            ["implicit_genborn_params", "cmaptypes", "cmaptypes"],
        )
        self.assertEqual(
            [section.name for section in reparsed.molecules["CMAP"][0].raw_sections],
            ["cmap"],
        )

    def test_manual_backed_optional_and_variable_sections_parse(self) -> None:
        topology_text = """
[ defaults ]
1 2

[ atomtypes ]
; name  bonded  at.num  mass   charge ptype sigma epsilon
C      12.011  0.0     A      0.34   0.10
H      HC      1.008   0.0    A      0.25   0.05
O      OT      8       15.999 -0.5   A      0.30   0.20
N      7       14.007  -0.3   A      0.32   0.15

[ constrainttypes ]
C H 1 0.109

[ moleculetype ]
MOL 3

[ atoms ]
1 C 1 MOL C1 1 0.0 12.011
2 H 1 MOL H1 2 0.0 1.008
3 O 1 MOL O1 3 -0.5 15.999
4 N 1 MOL N1 4 -0.3 14.007

[ constraints ]
1 2 1

[ angles ]
1 2 3 9 0.25 1000
2 3 4 10 120.0 50.0

[ virtual_sitesn ]
4 1 1 2 3

[ system ]
ManualCoverage

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "manual-coverage.top"
            path.write_text(topology_text)
            top = Topology(path)

        self.assertEqual(top.defaults.gen_pairs, "no")
        self.assertEqual(top.defaults.fudgeLJ, 1.0)
        self.assertEqual(top.defaults.fudgeQQ, 1.0)
        self.assertIsNone(top.defaults.n)

        atomtypes = {atomtype.name: atomtype for atomtype in top.atomtypes}
        self.assertIsNone(atomtypes["C"].bonded_type)
        self.assertIsNone(atomtypes["C"].atnum)
        self.assertEqual(atomtypes["H"].bonded_type, "HC")
        self.assertIsNone(atomtypes["H"].atnum)
        self.assertEqual(atomtypes["O"].bonded_type, "OT")
        self.assertEqual(atomtypes["O"].atnum, 8)
        self.assertIsNone(atomtypes["N"].bonded_type)
        self.assertEqual(atomtypes["N"].atnum, 7)

        molecule = top.molecules["MOL"][0]
        self.assertEqual(len(molecule.constraints), 1)
        self.assertEqual(molecule.constraints[0].params["b0"], 0.109)
        self.assertEqual(molecule.angles[0].func, 9)
        self.assertEqual(molecule.angles[0].params["a0"], 0.25)
        self.assertEqual(molecule.angles[1].func, 10)
        self.assertEqual(molecule.angles[1].params["th0"], 120.0)
        self.assertEqual(molecule.virtual_sitesn[0].params["from"], "1 2 3")

        with TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "manual-coverage-output.top"
            top.write(output)
            reparsed = Topology(output)

        self.assertEqual(len(reparsed.constrainttypes), 0)
        self.assertEqual(
            reparsed.molecules["MOL"][0].constraints[0].params["b0"],
            0.109,
        )

    def test_written_topology_is_flattened_and_filters_global_parameters(self) -> None:
        topology_text = """
[ defaults ]
1 2 no

[ atomtypes ]
A 12.0 0.0 A 0.30 0.10
B 14.0 0.0 A 0.32 0.12
UNUSED 16.0 0.0 A 0.34 0.14
RAW 18.0 0.0 A 0.36 0.16

[ nonbond_params ]
A B 1 0.31 0.11
A UNUSED 1 0.33 0.13

[ bondtypes ]
A B 1 0.10 1000

[ pairtypes ]
A B 1 0.20 0.30

[ angletypes ]
A B A 1 109.0 200

[ dihedraltypes ]
A B A B 9 180.0 4.0 3

[ constrainttypes ]
A A 1 0.20

[ opaque_global ]
RAW preserved

[ moleculetype ]
MOL 3

[ atoms ]
1 A 1 MOL A1 1 0.0 12.0
2 B 1 MOL B1 2 0.0 14.0
3 A 1 MOL A2 3 0.0 12.0
4 B 1 MOL B2 4 0.0 14.0

[ bonds ]
1 2 1

[ pairs ]
1 2 1

[ angles ]
1 2 3 1

[ dihedrals ]
1 2 3 4 9

[ constraints ]
1 3 1

[ system ]
Flattened

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "global-types.top"
            output = Path(tmpdir) / "flattened.top"
            path.write_text(topology_text)

            top = Topology(path)
            top.write(output)
            written = output.read_text()
            reparsed = Topology(output)

        self.assertEqual(
            {atomtype.name for atomtype in reparsed.atomtypes},
            {"A", "B", "RAW"},
        )
        self.assertEqual(len(reparsed.nonbond_params), 1)
        self.assertEqual(reparsed.nonbond_params[0].ai.name, "A")
        self.assertEqual(reparsed.nonbond_params[0].aj.name, "B")

        for section in (
            "bondtypes",
            "pairtypes",
            "angletypes",
            "dihedraltypes",
            "constrainttypes",
        ):
            self.assertNotIn(f"[ {section} ]", written)

        molecule = reparsed.molecules["MOL"][0]
        self.assertEqual(molecule.bonds[0].params["b0"], 0.10)
        self.assertEqual(molecule.pairs[0].params["epsilon"], 0.30)
        self.assertEqual(molecule.angles[0].params["th0"], 109.0)
        self.assertEqual(molecule.dihedrals[0].params["mult"], 3)
        self.assertEqual(molecule.constraints[0].params["b0"], 0.20)

    def test_filtered_output_preserves_conditional_section_positions(self) -> None:
        topology_text = """
[ defaults ]
1 2 no

[ atomtypes ]
#ifdef MODIFIED
A 12.0 0.0 A 0.30 0.10
B 14.0 0.0 A 0.32 0.12
UNUSED 16.0 0.0 A 0.34 0.14
#else
A 12.0 0.0 A 0.31 0.11
B 14.0 0.0 A 0.33 0.13
UNUSED 16.0 0.0 A 0.35 0.15
#endif

[ nonbond_params ]
#ifdef MODIFIED
A B 1 0.31 0.11
A UNUSED 1 0.33 0.13
#else
A B 1 0.32 0.12
A UNUSED 1 0.34 0.14
#endif

[ moleculetype ]
MOL 3

[ atoms ]
1 A 1 MOL A1 1 0.0 12.0
2 B 1 MOL B1 2 0.0 14.0

[ system ]
ConditionalFiltering

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conditional-filtering.top"
            output = Path(tmpdir) / "flattened.top"
            path.write_text(topology_text)

            top = Topology(path)
            top.write(output)
            written = output.read_text()
            reparsed = Topology(output)

        self.assertNotIn("UNUSED", written)
        self.assertEqual(written.count("#ifdef MODIFIED\n"), 2)
        self.assertEqual(written.count("#else\n"), 2)
        self.assertEqual(written.count("#endif\n"), 2)

        atomtypes_start = written.index("#ifdef MODIFIED\n")
        atomtypes_else = written.index("#else\n", atomtypes_start)
        atomtypes_end = written.index("#endif\n", atomtypes_else)
        nonbond_start = written.index("#ifdef MODIFIED\n", atomtypes_end)
        nonbond_else = written.index("#else\n", nonbond_start)
        nonbond_end = written.index("#endif\n", nonbond_else)

        self.assertLess(atomtypes_start, written.index("[ atomtypes ]"))
        self.assertLess(written.index("[ atomtypes ]"), atomtypes_else)
        self.assertLess(atomtypes_else, atomtypes_end)
        self.assertLess(atomtypes_end, nonbond_start)
        self.assertLess(nonbond_start, written.index("[ nonbond_params ]"))
        self.assertLess(written.index("[ nonbond_params ]"), nonbond_else)
        self.assertLess(nonbond_else, nonbond_end)

        expected_states = [
            ("ifdef MODIFIED",),
            ("ifdef MODIFIED",),
            ("else ifdef MODIFIED",),
            ("else ifdef MODIFIED",),
        ]
        self.assertEqual(
            [atomtype.ifdef_state for atomtype in reparsed.atomtypes],
            expected_states,
        )
        self.assertEqual(
            [param.ifdef_state for param in reparsed.nonbond_params],
            [("ifdef MODIFIED",), ("else ifdef MODIFIED",)],
        )

    def test_filtered_output_reconstructs_empty_conditional_branches(self) -> None:
        topology_text = """
[ defaults ]
1 2 no

[ atomtypes ]
A 12.0 0.0 A 0.30 0.10
B 14.0 0.0 A 0.32 0.12
UNUSED 16.0 0.0 A 0.34 0.14

[ nonbond_params ]
#ifdef OUTER
#ifdef OMITTED
A UNUSED 1 0.33 0.13
#else
A B 1 0.31 0.11
#endif
#endif

[ moleculetype ]
MOL 3

[ atoms ]
1 A 1 MOL A1 1 0.0 12.0
2 B 1 MOL B1 2 0.0 14.0

[ system ]
EmptyConditionalBranch

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conditional-filtering.top"
            output = Path(tmpdir) / "flattened.top"
            path.write_text(topology_text)

            top = Topology(path)
            top.write(output)
            written = output.read_text()
            reparsed = Topology(output)

        self.assertNotIn("UNUSED", written)
        self.assertIn(
            "#ifdef OUTER\n"
            "#ifdef OMITTED\n"
            "#else\n"
            "[ nonbond_params ]",
            written,
        )
        self.assertIn("[ nonbond_params ]", written)
        self.assertIn("#endif\n#endif\n", written)
        self.assertEqual(
            reparsed.nonbond_params[0].ifdef_state,
            ("ifdef OUTER", "else ifdef OMITTED"),
        )

    def test_conditional_defines_roundtrip_in_their_original_branches(self) -> None:
        topology_text = """
[ defaults ]
1 2 no

#ifdef HEAVY
#define PARTICLE_MASS 4.0
#else
#define PARTICLE_MASS 1.0
#endif

[ atomtypes ]
A PARTICLE_MASS 0.0 A 0.30 0.10

[ moleculetype ]
MOL 3

[ atoms ]
1 A 1 MOL A1 1 0.0 PARTICLE_MASS

[ system ]
ConditionalDefines

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conditional-defines.top"
            output = Path(tmpdir) / "flattened.top"
            path.write_text(topology_text)

            top = Topology(path)
            top.write(output)
            written = output.read_text()
            reparsed = Topology(output)

        self.assertIn(
            "#ifdef HEAVY\n"
            "#define PARTICLE_MASS 4.0\n"
            "#else\n"
            "#define PARTICLE_MASS 1.0\n"
            "#endif\n",
            written,
        )
        self.assertEqual(
            [define.ifdef_state for define in reparsed.defines],
            [("ifdef HEAVY",), ("else ifdef HEAVY",)],
        )

    def test_marker_define_and_ifndef_roundtrip(self) -> None:
        topology_text = """
#define POSRES
#define RESTRAINT_FC 1000 1000 1000

[ defaults ]
1 2

[ atomtypes ]
C 12.011 0.0 A 0.34 0.10

[ moleculetype ]
MOL 3

[ atoms ]
#ifndef FLEXIBLE
1 C 1 MOL C1 1 0.0 12.011
#endif

[ system ]
Conditional

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conditional.top"
            output = Path(tmpdir) / "conditional-output.top"
            path.write_text(topology_text)

            top = Topology(path)
            defines = {define.directive: define.argument for define in top.defines}
            atom = top.molecules["MOL"][0].atoms[0]

            self.assertIsNone(defines["POSRES"])
            self.assertEqual(defines["RESTRAINT_FC"], "1000 1000 1000")
            self.assertEqual(atom.ifdef_state, ("ifndef FLEXIBLE",))

            top.write(output)
            written = output.read_text()
            self.assertIn("#define POSRES\n", written)
            self.assertIn("#define RESTRAINT_FC 1000 1000 1000\n", written)
            self.assertLess(
                written.index("#define POSRES\n"),
                written.index("#define RESTRAINT_FC 1000 1000 1000\n"),
            )
            self.assertIn("#ifndef FLEXIBLE\n", written)

            reparsed = Topology(output)
            reparsed_atom = reparsed.molecules["MOL"][0].atoms[0]
            self.assertEqual(reparsed_atom.ifdef_state, ("ifndef FLEXIBLE",))

    def test_free_energy_b_state_parameters_roundtrip(self) -> None:
        topology_text = """
[ defaults ]
1 2

[ atomtypes ]
A 12.0 0.0 A 0.34 0.10
B 14.0 0.0 A 0.36 0.12

[ moleculetype ]
MOL 3

[ atoms ]
1 A 1 MOL A1 1 0.0 B -0.2 14.0
2 A 1 MOL A2 1 0.0
3 A 1 MOL A3 1 0.0 12.0

[ bonds ]
1 2 1 0.1 345000 0.2 300000

[ angles ]
1 2 3 1 109.47 383 120.0 400

[ system ]
FreeEnergy

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "free-energy.top"
            output = Path(tmpdir) / "free-energy-output.top"
            path.write_text(topology_text)

            top = Topology(path)
            molecule = top.molecules["MOL"][0]

            self.assertIsNone(molecule.atoms[0].mass)
            self.assertIsNone(molecule.atoms[1].mass)
            self.assertEqual(molecule.atoms[0].type_b.name, "B")
            self.assertEqual(molecule.atoms[0].charge_b, -0.2)
            self.assertEqual(molecule.atoms[0].mass_b, 14.0)
            self.assertEqual(molecule.bonds[0].params["b0_b"], 0.2)
            self.assertEqual(molecule.bonds[0].params["kb_b"], 300000.0)
            self.assertEqual(molecule.angles[0].params["th0_b"], 120.0)
            self.assertEqual(molecule.angles[0].params["kth_b"], 400.0)

            top.write(output)
            reparsed = Topology(output)

        reparsed_molecule = reparsed.molecules["MOL"][0]
        self.assertEqual(reparsed_molecule.atoms[0].type_b.name, "B")
        self.assertEqual(reparsed_molecule.bonds[0].params["kb_b"], 300000.0)
        self.assertEqual(reparsed_molecule.angles[0].params["kth_b"], 400.0)

    def test_nested_conditional_include_roundtrip(self) -> None:
        topology_text = """
[ defaults ]
1 2

[ atomtypes ]
C 12.011 0.0 A 0.34 0.10

#ifndef FLEXIBLE
#include "molecule.itp"
#endif

[ system ]
Included

[ molecules ]
MOL 1
""".strip()
        molecule_text = """
[ moleculetype ]
MOL 3

[ atoms ]
#ifdef HEAVY
1 C 1 MOL C1 1 0.0 12.011
#else
1 C 1 MOL C1 1 0.0 12.011
#endif
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "conditional-include.top"
            include = Path(tmpdir) / "molecule.itp"
            output = Path(tmpdir) / "conditional-include-output.top"
            path.write_text(topology_text)
            include.write_text(molecule_text)

            top = Topology(path)
            molecule = top.molecules["MOL"][0]
            self.assertEqual(molecule.ifdef_state, ("ifndef FLEXIBLE",))
            self.assertEqual(
                [atom.ifdef_state for atom in molecule.atoms],
                [
                    ("ifndef FLEXIBLE", "ifdef HEAVY"),
                    ("ifndef FLEXIBLE", "else ifdef HEAVY"),
                ],
            )

            top.write(output)
            reparsed = Topology(output)
            reparsed_molecule = reparsed.molecules["MOL"][0]
            self.assertEqual(reparsed_molecule.ifdef_state, ("ifndef FLEXIBLE",))
            self.assertEqual(
                [atom.ifdef_state for atom in reparsed_molecule.atoms],
                [
                    ("ifndef FLEXIBLE", "ifdef HEAVY"),
                    ("ifndef FLEXIBLE", "else ifdef HEAVY"),
                ],
            )

    def test_missing_include_is_allowed_only_inside_a_conditional(self) -> None:
        unconditional_text = """
[ defaults ]
1 2

#include "missing.itp"
""".strip()
        conditional_text = """
[ defaults ]
1 2

#ifdef OPTIONAL
#include "missing.itp"
#endif

[ atomtypes ]
C 12.011 0.0 A 0.34 0.10

[ moleculetype ]
MOL 3

[ atoms ]
1 C 1 MOL C1 1 0.0 12.011

[ system ]
MissingOptionalInclude

[ molecules ]
MOL 1
""".strip()

        with TemporaryDirectory() as tmpdir:
            unconditional = Path(tmpdir) / "unconditional.top"
            conditional = Path(tmpdir) / "conditional.top"
            unconditional.write_text(unconditional_text)
            conditional.write_text(conditional_text)

            with self.assertRaises(FileNotFoundError):
                Topology(unconditional)

            top = Topology(conditional)

        self.assertEqual(len(top.atoms), 1)

    def test_molecule_type_names_are_case_insensitive(self) -> None:
        topology_text = """
[ defaults ]
1 2

[ atomtypes ]
C 12.011 0.0 A 0.34 0.10

[ moleculetype ]
MOL 3

[ atoms ]
1 C 1 MOL C1 1 0.0 12.011

[ system ]
MoleculeNameCase

[ molecules ]
mol 1
""".strip()
        duplicate_text = """
[ defaults ]
1 2

[ moleculetype ]
MOL 3

[ moleculetype ]
mol 3
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "case-insensitive.top"
            duplicate = Path(tmpdir) / "duplicate.top"
            path.write_text(topology_text)
            duplicate.write_text(duplicate_text)

            top = Topology(path)

            with self.assertRaisesRegex(ValueError, "Duplicate moleculetype"):
                Topology(duplicate)

        self.assertEqual(list(top.molecules), ["MOL"])

    def test_molecule_fragment_include_inherits_context(self) -> None:
        topology_text = """
[ defaults ]
1 2

[ atomtypes ]
C 12.011 0.0 A 0.34 0.10

[ moleculetype ]
MOL 3

[ atoms ]
1 C 1 MOL C1 1 0.0 12.011

#ifdef POSRES
#include "posre.itp"
#endif

[ system ]
IncludedRestraints

[ molecules ]
MOL 1
""".strip()
        restraint_text = """
[ position_restraints ]
1 1 1000 1000 1000
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "topol.top"
            include = Path(tmpdir) / "posre.itp"
            path.write_text(topology_text)
            include.write_text(restraint_text)

            top = Topology(path)

        restraints = top.molecules["MOL"][0].position_restraints
        self.assertEqual(len(restraints), 1)
        self.assertEqual(restraints[0].params["kx"], 1000.0)
        self.assertEqual(restraints[0].ifdef_state, ("ifdef POSRES",))

    def test_unclosed_conditional_is_rejected(self) -> None:
        topology_text = """
[ defaults ]
1 2

#ifndef FLEXIBLE
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "unclosed.top"
            path.write_text(topology_text)

            with self.assertRaisesRegex(ValueError, "Unclosed conditional block"):
                Topology(path)

    def test_unsupported_preprocessor_directive_is_rejected(self) -> None:
        topology_text = """
[ defaults ]
1 2

#undef POSRES
""".strip()

        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "unsupported-directive.top"
            path.write_text(topology_text)

            with self.assertRaisesRegex(
                NotImplementedError,
                "Unsupported preprocessor directive",
            ):
                Topology(path)


if __name__ == "__main__":
    unittest.main()
