"""GROMACS interaction schemas keyed by topology directive and function type."""

from .schema import ParamSpec, InteractionSpec


BONDS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="harmonic bond",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="G96 bond",
            allow_b_state=True,
        ),
        3: ParamSpec(
            names=("b0", "D", "beta"),
            parsers=(float, float, float),
            desc="Morse bond",
            allow_b_state=True,
        ),
        4: ParamSpec(
            names=("b0", "C2", "C3"),
            parsers=(float, float, float),
            desc="cubic bond",
        ),
        5: ParamSpec(names=(), parsers=(), desc="connection (no parameters)"),
        6: ParamSpec(
            names=("b0", "kb"),
            parsers=(float, float),
            desc="harmonic potential",
            allow_b_state=True,
        ),
        7: ParamSpec(names=("bm", "kb"), parsers=(float, float), desc="FENE bond"),
        8: ParamSpec(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated bond",
            allow_b_state=True,
        ),
        9: ParamSpec(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated bond (variant)",
            allow_b_state=True,
        ),
        10: ParamSpec(
            names=("low", "up1", "up2", "kdr"),
            parsers=(float, float, float, float),
            desc="restraint",
            allow_b_state=True,
        ),
    },
)

PAIRS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("sigma", "epsilon"),
            parsers=(float, float),
            desc="Lennard-Jones pair",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("fudgeQQ", "q1", "q2", "sigma", "epsilon"),
            parsers=(float, float, float, float, float),
            desc="Coulomb + LJ pair",
        ),
    },
)

PAIRS_NB = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("qi", "qj", "V", "W"),
            parsers=(float, float, float, float),
            desc="non-bonded pair",
        ),
    },
)

ANGLES = InteractionSpec(
    n_atoms=3,
    funcs={
        1: ParamSpec(
            names=("th0", "kth"),
            parsers=(float, float),
            desc="angle",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("th0", "kth"),
            parsers=(float, float),
            desc="G96 angle",
            allow_b_state=True,
        ),
        3: ParamSpec(
            names=("r1e", "r2e", "krr"),
            parsers=(float, float, float),
            desc="cross bond-bond",
        ),
        4: ParamSpec(
            names=("r1e", "r2e", "r3e", "krth"),
            parsers=(float, float, float, float),
            desc="cross bond-angle",
        ),
        5: ParamSpec(
            names=("th0", "kth", "r13", "kub"),
            parsers=(float, float, float, float),
            desc="Urey-Bradley",
            allow_b_state=True,
        ),
        6: ParamSpec(
            names=("th0", "C0", "C1", "C2", "C3", "C4"),
            parsers=(float, float, float, float, float, float),
            desc="quartic angle",
        ),
        8: ParamSpec(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated angle",
            allow_b_state=True,
        ),
        9: ParamSpec(
            names=("a0", "klin"),
            parsers=(float, float),
            desc="linear angle",
            allow_b_state=True,
        ),
        10: ParamSpec(
            names=("th0", "kth"),
            parsers=(float, float),
            desc="restricted bending potential",
        ),
    },
)

DIHEDRALS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: ParamSpec(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="proper dihedral",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("xi0", "kxi"),
            parsers=(float, float),
            desc="improper dihedral",
            allow_b_state=True,
        ),
        3: ParamSpec(
            names=("C0", "C1", "C2", "C3", "C4", "C5"),
            parsers=(float, float, float, float, float, float),
            desc="RB dihedral",
            allow_b_state=True,
        ),
        4: ParamSpec(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="periodic improper",
            allow_b_state=True,
        ),
        5: ParamSpec(
            names=("C1", "C2", "C3", "C4", "C5"),
            parsers=(float, float, float, float, float),
            desc="Fourier dihedral",
            allow_b_state=True,
        ),
        8: ParamSpec(
            names=("table", "k"),
            parsers=(int, float),
            desc="tabulated dihedral",
            allow_b_state=True,
        ),
        9: ParamSpec(
            names=("phi_s", "kphi", "mult"),
            parsers=(float, float, int),
            desc="proper dihedral (multiple)",
            allow_b_state=True,
        ),
        10: ParamSpec(
            names=("phi0", "kphi"),
            parsers=(float, float),
            desc="restricted dihedral",
        ),
        11: ParamSpec(
            names=("kphi", "a0", "a1", "a2", "a3", "a4"),
            parsers=(float, float, float, float, float, float),
            desc="combined bending-torsion",
        ),
    },
)

CONSTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("b0",),
            parsers=(float,),
            desc="constraint",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("b0",),
            parsers=(float,),
            desc="constraint",
            allow_b_state=True,
        ),
    },
)

NONBOND_PARAMS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("sigma", "epsilon"),
            parsers=(float, float),
            desc="Lennard-Jones parameters",
        ),
        2: ParamSpec(
            names=("a", "b", "c6"),
            parsers=(float, float, float),
            desc="Buckingham parameters",
        ),
    },
)

SETTLES = InteractionSpec(
    n_atoms=1,
    funcs={
        1: ParamSpec(
            names=("doh", "dhh"),
            parsers=(float, float),
            desc="SETTLES parameters",
        ),
    },
)

VIRTUAL_SITES1 = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(names=("a",), parsers=(float,), desc="1-body virtual site"),
    },
)

VIRTUAL_SITES2 = InteractionSpec(
    n_atoms=3,
    funcs={
        1: ParamSpec(names=("a",), parsers=(float,), desc="2-body virtual site"),
        2: ParamSpec(
            names=("d",),
            parsers=(float,),
            desc="2-body virtual site (fd)",
        ),
    },
)

VIRTUAL_SITES3 = InteractionSpec(
    n_atoms=4,
    funcs={
        1: ParamSpec(
            names=("a", "b"),
            parsers=(float, float),
            desc="3-body virtual site",
        ),
        2: ParamSpec(
            names=("a", "d"),
            parsers=(float, float),
            desc="3-body virtual site (fd)",
        ),
        3: ParamSpec(
            names=("th", "d"),
            parsers=(float, float),
            desc="3-body virtual site (fad)",
        ),
        4: ParamSpec(
            names=("a", "b", "c"),
            parsers=(float, float, float),
            desc="3-body virtual site (out)",
        ),
    },
)

VIRTUAL_SITES4 = InteractionSpec(
    n_atoms=5,
    funcs={
        2: ParamSpec(
            names=("a", "b", "c"),
            parsers=(float, float, float),
            desc="4-body virtual site (fdn)",
        ),
    },
)

VIRTUAL_SITESN = InteractionSpec(
    n_atoms=1,
    funcs={
        1: ParamSpec(
            names=(),
            parsers=(),
            desc="N-body virtual site (COG)",
            rest_name="from",
        ),
        2: ParamSpec(
            names=(),
            parsers=(),
            desc="N-body virtual site (COM)",
            rest_name="from",
        ),
        3: ParamSpec(
            names=(),
            parsers=(),
            desc="N-body virtual site (COW)",
            rest_name="from",
        ),
    },
)

POSITION_RESTRAINTS = InteractionSpec(
    n_atoms=1,
    funcs={
        1: ParamSpec(
            names=("kx", "ky", "kz"),
            parsers=(float, float, float),
            desc="position restraints",
            allow_b_state=True,
        ),
        2: ParamSpec(
            names=("g", "r", "k"),
            parsers=(float, float, float),
            desc="position restraints variant",
        ),
    },
)

DISTANCE_RESTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("type", "label", "low", "up1", "up2", "weight"),
            parsers=(int, int, float, float, float, float),
            desc="distance restraints",
        ),
    },
)

DIHEDRAL_RESTRAINTS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: ParamSpec(
            names=("phi0", "dphi", "kdihr"),
            parsers=(float, float, int),
            desc="dihedral restraints",
            allow_b_state=True,
        ),
    },
)

ORIENTATION_RESTRAINTS = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("exp", "label", "alpha", "c", "obs", "weight"),
            parsers=(int, int, float, float, float, float),
            desc="orientation restraints",
        ),
    },
)

ANGLE_RESTRAINTS = InteractionSpec(
    n_atoms=4,
    funcs={
        1: ParamSpec(
            names=("theta0", "kc", "mult"),
            parsers=(float, float, int),
            desc="angle restraints",
            allow_b_state=True,
        ),
    },
)

ANGLE_RESTRAINTS_Z = InteractionSpec(
    n_atoms=2,
    funcs={
        1: ParamSpec(
            names=("theta0", "kc", "mult"),
            parsers=(float, float, int),
            desc="angle restraints z",
            allow_b_state=True,
        ),
    },
)
