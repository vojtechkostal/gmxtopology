"""GROMACS interaction directives and the parameters of each function type.

``SPECS`` maps a directive name to the number of atoms it takes and, for each
function type, the names of its parameters. A ``:int`` suffix marks integer
parameters (all others are floats). Any function type may be followed by a
second (B-state) set of the same parameters, stored with a ``_b`` suffix.
"""

from dataclasses import dataclass

from .preprocess import TopologyError


@dataclass(frozen=True)
class Func:
    """Parameters of one function type of an interaction directive."""

    names: tuple[str, ...]
    ints: frozenset[str] = frozenset()
    rest: str | None = None  # name under which a variable-length tail is stored


@dataclass(frozen=True)
class Spec:
    n_atoms: int
    funcs: dict[int, Func]


def _f(names: str = "", *, rest: str | None = None) -> Func:
    tokens = [token.partition(":") for token in names.split()]
    return Func(
        names=tuple(name for name, _, _ in tokens),
        ints=frozenset(name for name, _, kind in tokens if kind == "int"),
        rest=rest,
    )


SPECS: dict[str, Spec] = {
    "bonds": Spec(
        2,
        {
            1: _f("b0 kb"),
            2: _f("b0 kb"),
            3: _f("b0 D beta"),
            4: _f("b0 C2 C3"),
            5: _f(),
            6: _f("b0 kb"),
            7: _f("bm kb"),
            8: _f("table:int k"),
            9: _f("table:int k"),
            10: _f("low up1 up2 kdr"),
        },
    ),
    "pairs": Spec(
        2,
        {
            1: _f("sigma epsilon"),
            2: _f("fudgeQQ qi qj sigma epsilon"),
        },
    ),
    "pairs_nb": Spec(2, {1: _f("qi qj V W")}),
    "angles": Spec(
        3,
        {
            1: _f("th0 kth"),
            2: _f("th0 kth"),
            3: _f("r1e r2e krr"),
            4: _f("r1e r2e r3e krth"),
            5: _f("th0 kth r13 kub"),
            6: _f("th0 C0 C1 C2 C3 C4"),
            8: _f("table:int k"),
            9: _f("a0 klin"),
            10: _f("th0 kth"),
        },
    ),
    "dihedrals": Spec(
        4,
        {
            1: _f("phi_s kphi mult:int"),
            2: _f("xi0 kxi"),
            3: _f("C0 C1 C2 C3 C4 C5"),
            4: _f("phi_s kphi mult:int"),
            5: _f("C1 C2 C3 C4"),
            8: _f("table:int k"),
            9: _f("phi_s kphi mult:int"),
            10: _f("phi0 kphi"),
            11: _f("kphi a0 a1 a2 a3 a4"),
        },
    ),
    "constraints": Spec(
        2,
        {
            1: _f("b0"),
            2: _f("b0"),
        },
    ),
    "settles": Spec(1, {1: _f("doh dhh")}),
    "virtual_sites1": Spec(2, {1: _f()}),
    "virtual_sites2": Spec(3, {1: _f("a"), 2: _f("d")}),
    "virtual_sites3": Spec(
        4,
        {
            1: _f("a b"),
            2: _f("a d"),
            3: _f("th d"),
            4: _f("a b c"),
        },
    ),
    "virtual_sites4": Spec(5, {2: _f("a b c")}),
    "virtual_sitesn": Spec(
        1,
        {
            1: _f(rest="from"),
            2: _f(rest="from"),
            3: _f(rest="from"),
        },
    ),
    "position_restraints": Spec(
        1,
        {
            1: _f("kx ky kz"),
            2: _f("g r k"),
        },
    ),
    "distance_restraints": Spec(
        2,
        {
            1: _f("type:int label:int low up1 up2 weight"),
        },
    ),
    "dihedral_restraints": Spec(4, {1: _f("phi0 dphi kdihr")}),
    "orientation_restraints": Spec(
        2,
        {
            1: _f("exp:int label:int alpha c obs weight"),
        },
    ),
    "angle_restraints": Spec(4, {1: _f("theta0 kc mult:int")}),
    "angle_restraints_z": Spec(2, {1: _f("theta0 kc mult:int")}),
    "cmap": Spec(5, {1: _f()}),
    "nonbond_params": Spec(
        2,
        {
            1: _f("sigma epsilon"),
            2: _f("a b c6"),
        },
    ),
}

# Parameter tables (``[ bondtypes ]`` ...) share the parameters of the molecular
# directive they provide defaults for.
TYPE_SECTIONS: dict[str, str] = {
    "bondtypes": "bonds",
    "pairtypes": "pairs",
    "angletypes": "angles",
    "dihedraltypes": "dihedrals",
    "constrainttypes": "constraints",
    "nonbond_params": "nonbond_params",
}

# Directives that take a default function type of 1 when only atoms are given.
DEFAULT_FUNC_SECTIONS = frozenset(
    {"bonds", "pairs", "angles", "dihedrals", "constraints"}
)

# Directives whose parameters may come from a parameter table.
LOOKUP_SECTIONS = frozenset({"bonds", "pairs", "angles", "dihedrals", "constraints"})

ALIASES = {
    "dummies1": "virtual_sites1",
    "dummies2": "virtual_sites2",
    "dummies3": "virtual_sites3",
    "dummies4": "virtual_sites4",
    "dummiesn": "virtual_sitesn",
}

VSITE_SECTIONS = tuple(name for name in SPECS if name.startswith("virtual_sites"))
MOLECULE_SECTIONS = tuple(name for name in SPECS if name not in {"nonbond_params"})
INTERMOLECULAR_SECTIONS = ("bonds", "angles", "dihedrals")


def parse_value(token: str, kind: type, macros, where: str) -> float | int | str:
    """Convert a token to ``kind``, or keep it as a macro name if it is one."""
    try:
        if kind is int:
            number = float(token)
            if number != int(number):
                raise ValueError
            return int(number)
        return float(token)
    except ValueError:
        if token.lstrip("-") in macros:
            return token
        raise TopologyError(
            f"Invalid {kind.__name__} value '{token}' at {where}."
        ) from None


def parse_params(
    section: str, func: int, tokens: list[str], macros, where: str
) -> dict[str, float | int | str]:
    """Parse the parameter tokens of one interaction line."""
    spec = SPECS[section]
    if func not in spec.funcs:
        raise TopologyError(
            f"Function type {func} is not supported for [ {section} ] at {where}; "
            f"supported types are {', '.join(map(str, spec.funcs))}."
        )
    fn = spec.funcs[func]
    n = len(fn.names)
    if len(tokens) < n:
        raise TopologyError(
            f"Expected {n} parameter(s) for [ {section} ] function {func} at "
            f"{where}, got {len(tokens)}."
        )

    groups = [("", tokens[:n])]
    rest = tokens[n:]
    if n and len(rest) == n:
        groups.append(("_b", rest))
        rest = []

    params: dict[str, float | int | str] = {}
    for suffix, values in groups:
        for name, token in zip(fn.names, values):
            kind = int if name in fn.ints else float
            params[name + suffix] = parse_value(token, kind, macros, where)

    if fn.rest is not None:
        if not rest:
            raise TopologyError(
                f"Expected atoms after the function type of [ {section} ] at {where}."
            )
        params[fn.rest] = " ".join(rest)
    elif rest:
        raise TopologyError(
            f"Too many parameters for [ {section} ] function {func} at {where}: "
            f"{' '.join(rest)}."
        )
    return params
