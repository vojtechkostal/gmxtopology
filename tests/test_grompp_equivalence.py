"""A written topology must describe the same system as the one it was read from.

Both topologies are run through ``gmx grompp`` and the interactions of the
resulting run inputs are compared with ``gmx dump``. The tests are skipped when
GROMACS is not installed.
"""

import random
import re
import shutil
import subprocess
import warnings
from pathlib import Path

import pytest

from gmxtopology import Topology

GMX = shutil.which("gmx")
pytestmark = pytest.mark.skipif(GMX is None, reason="GROMACS (gmx) is not installed")

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures"

MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
pbc = xyz
rvdw = 0.5
rcoulomb = 0.5
define = {define}
include = -I{include}
"""

PEPTIDE = {
    "ALA": [("CB", 0, -1.5, 0.5)],
    "GLU": [
        ("CB", 0, -1.5, 0.5),
        ("CG", 0, -2.9, 1.0),
        ("CD", 0, -4.2, 1.5),
        ("OE1", 0, -5.2, 0.9),
        ("OE2", 0, -4.3, 2.7),
    ],
    "SER": [("CB", 0, -1.5, 0.5), ("OG", 0, -2.7, 1.0)],
}


def gmx(*args: str, cwd: Path) -> str:
    result = subprocess.run([GMX, *args], cwd=cwd, capture_output=True, text=True)
    return result.stdout + result.stderr


def write_gro(top: Topology, path: Path) -> None:
    """Random coordinates in a big box; grompp only needs matching atom counts."""
    rng = random.Random(1)
    atoms = top.atoms
    rows = [
        f"{a.resnr % 100000:5d}{a.resname[:5]:<5s}{a.name[:5]:>5s}{(i + 1) % 100000:5d}"
        + "".join(f"{rng.uniform(1, 49):8.3f}" for _ in range(3))
        for i, a in enumerate(atoms)
    ]
    path.write_text("\n".join(["t", str(len(atoms)), *rows, "50 50 50", ""]))


def run_grompp(
    top_path: Path, gro: Path, work: Path, name: str, define: str, include: Path
):
    (work / "p.mdp").write_text(MDP.format(define=define, include=include))
    out = gmx(
        "grompp",
        "-f",
        "p.mdp",
        "-c",
        str(gro),
        "-r",
        str(gro),
        "-p",
        str(top_path),
        "-o",
        f"{name}.tpr",
        "-maxwarn",
        "50",
        cwd=work,
    )
    assert (work / f"{name}.tpr").exists(), (
        f"grompp failed for {top_path}:\n{out[-2000:]}"
    )
    return gmx("dump", "-s", f"{name}.tpr", cwd=work)


def summarize(dump: str) -> dict:
    """Interactions with their resolved parameters, independent of table order."""
    text = dump[dump.index("topology:") :]
    functype = {int(i): body for i, body in re.findall(r"functype\[(\d+)\]=(.*)", text)}
    grids = [
        re.split(r"\n\s*V\s+dVdx|\n   moltype", p)[0].strip()
        for p in re.split(r"grid\[\s*\d+\]=\{", text)[1:]
    ]
    for i, body in functype.items():
        m = re.match(r"CMAP, cmapA=(\d+)", body)
        if m:
            functype[i] = "CMAP " + grids[int(m.group(1))]
    type_names = {
        int(i): name for i, name in re.findall(r"type\[(\d+)\]=\{name=\"(.*?)\"", text)
    }
    summary = {"lj": {}, "molecules": []}
    atnr = int(re.search(r"atnr=(\d+)", text).group(1))
    summary["lj"] = {
        (type_names[i], type_names[j]): functype[i * atnr + j]
        for i in range(atnr)
        for j in range(atnr)
    }
    for block in re.split(r"\n   moltype \(\d+\):", text)[1:]:
        block = block.split("\n   molblock")[0]
        atoms = [
            (type_names[int(t)], type_names[int(tb)], rest)
            for t, tb, rest in re.findall(
                r"atom\[\s*\d+\]=\{type=\s*(\d+), typeB=\s*(\d+), (.*)\}", block
            )
        ]
        ilist: dict[str, list] = {}
        section = None
        for line in block.splitlines():
            if m := re.match(r"      (\S.*):$", line):
                section = m.group(1)
            if m := re.match(r"\s+\d+ type=(\d+) \(.*?\)\s+(.*)", line):
                ilist.setdefault(section, []).append(
                    (functype[int(m.group(1))], m.group(2))
                )
        summary["molecules"].append(
            {
                "atoms": atoms,
                "excls": re.findall(r"excls\[\d+\]\[num=\d+\]=\{(.*?)\}", block),
                "ilist": {k: sorted(v) for k, v in ilist.items()},
            }
        )
    return summary


def assert_equivalent(top_path: Path, work: Path, defines: dict | None = None):
    defines = defines or {}
    define = " ".join(f"-D{k}" + (f"={v}" if v else "") for k, v in defines.items())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        top = Topology(top_path, defines=defines)
    gro = work / "p.gro"
    write_gro(top, gro)
    top.write(work / "flat.top")

    before = summarize(run_grompp(top_path, gro, work, "a", define, top_path.parent))
    after = summarize(run_grompp(work / "flat.top", gro, work, "b", define, work))
    assert after["molecules"] == before["molecules"]
    for pair, parameters in after["lj"].items():  # unused atom types are dropped
        assert before["lj"][pair] == parameters


@pytest.mark.parametrize(
    "relative",
    [
        "gromacs-v2026.2/topol-tip4p.top",
        "gromacs-v2026.2/topol-urea.top",
        "parmed-96ec61a/03.AlaGlu/topol.top",
        "parmed-96ec61a/12.DPPC/topol.top",
        "prosECCo75-e4831a4/topol-cmap.top",
        "prosECCo75-e4831a4/topol-popc.top",
    ],
)
def test_fixture_topologies_are_equivalent_under_grompp(tmp_path, relative):
    assert_equivalent(FIXTURES / relative, tmp_path)


@pytest.fixture(scope="module")
def peptide_pdb(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("peptide") / "pep.pdb"
    lines, n = [], 1
    residues = ["ALA", "GLU", "SER", "ALA", "ALA"]
    for i, res in enumerate(residues):
        x0 = i * 3.8
        atoms = [
            ("N", -1.3, 0, 0),
            ("CA", 0, 0, 0),
            ("C", 1.3, 0.5, 0),
            ("O", 1.4, 1.7, 0),
        ]
        atoms += [
            (nm, dx, dy, dz) for nm, dx, dy, dz in PEPTIDE.get(res, PEPTIDE["ALA"])
        ]
        if i == len(residues) - 1:
            atoms.append(("OXT", 2.4, -0.3, 0))
        for nm, dx, dy, dz in atoms:
            lines.append(
                f"ATOM  {n:5d} {nm:<4s} {res} A{i + 1:4d}    "
                f"{x0 + dx:8.3f}{dy:8.3f}{dz:8.3f}  1.00  0.00"
            )
            n += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


@pytest.mark.parametrize(
    "forcefield, defines",
    [
        ("amber99sb-ildn", {}),
        ("charmm27", {"POSRES": None}),  # CMAP and conditional position restraints
        ("oplsaa", {"FLEXIBLE": None}),  # bonded types differ from atom type names
        ("gromos54a7", {}),  # macros for bonded parameters
    ],
)
def test_pdb2gmx_topologies_are_equivalent_under_grompp(
    tmp_path, peptide_pdb, forcefield, defines
):
    water = "spc" if forcefield.startswith("gromos") else "tip3p"
    out = gmx(
        "pdb2gmx",
        "-f",
        str(peptide_pdb),
        "-ff",
        forcefield,
        "-water",
        water,
        "-ignh",
        "-o",
        "conf.gro",
        "-p",
        "topol.top",
        "-i",
        "posre.itp",
        cwd=tmp_path,
    )
    assert (tmp_path / "topol.top").exists(), out[-1500:]
    work = tmp_path / "work"
    work.mkdir()
    assert_equivalent(tmp_path / "topol.top", work, defines)
