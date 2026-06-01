from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence


@dataclass(frozen=True)
class ParamSpec:
    """Schema for the parameters of a single interaction function."""

    names: tuple[str, ...]
    parsers: tuple[type, ...]
    desc: str = ""
    rest_name: str | None = None
    rest_parser: Callable[[Sequence[str]], Any] | None = None
    allow_b_state: bool = False

    def __post_init__(self) -> None:
        if len(self.names) != len(self.parsers):
            raise ValueError(
                f"ParamSpec names and parsers must have the same length: "
                f"got {len(self.names)} names and {len(self.parsers)} parsers."
            )


@dataclass(frozen=True)
class InteractionSpec:
    """Schema for one interaction directive and its supported function types."""

    n_atoms: int
    funcs: Mapping[int, ParamSpec]

    def parse(self, func: int, tokens: Sequence[str], ctx: str = "") -> dict[str, Any]:
        """Parse tokens for one interaction line."""

        param_spec = self.funcs.get(func)
        if param_spec is None:
            allowed = ", ".join(map(str, self.funcs))
            raise NotImplementedError(
                f"Only function(s) {allowed} are supported in {ctx}: got {func}."
            )

        expected = len(param_spec.names)
        if len(tokens) < expected:
            raise ValueError(
                f"Expected {expected} parameter value(s) in {ctx} for function "
                f"{func}, got {len(tokens)}."
            )

        rest = tokens[expected:]
        parameter_sets = [("", tokens[:expected])]
        if param_spec.allow_b_state and len(rest) == expected:
            parameter_sets.append(("_b", rest))
            rest = ()

        params: dict[str, Any] = {}
        for suffix, values in parameter_sets:
            for name, parser, token in zip(
                param_spec.names,
                param_spec.parsers,
                values,
            ):
                try:
                    params[f"{name}{suffix}"] = parser(token)
                except ValueError:
                    params[f"{name}{suffix}"] = token

        if param_spec.rest_name is not None:
            if not rest:
                raise ValueError(
                    f"Expected additional parameter value(s) in {ctx} for function "
                    f"{func}."
                )
            params[param_spec.rest_name] = (
                param_spec.rest_parser(rest)
                if param_spec.rest_parser is not None
                else " ".join(rest)
            )
        elif rest:
            raise NotImplementedError(
                f"Extra parameter value(s) are not supported in {ctx} for function "
                f"{func}: {rest!r}."
            )

        return params
