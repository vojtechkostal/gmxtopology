"""Default parameters of molecular interactions from the parameter tables.

The rules follow grompp: bonded tables are matched on the *bonded* atom types
(either orientation), dihedrals also match wildcards (``X``) and the entry with
the most explicitly matched types wins, a ``[ dihedraltypes ]`` function-9
block is a run of adjacent lines with identical types, and a repeated entry of
any other kind replaces the earlier one.
"""

from collections import defaultdict
from typing import Iterable

from .model import TypeEntry

WILDCARD = "X"


def _canonical(types: tuple[str, ...]) -> tuple[str, ...]:
    return min(types, types[::-1])


def _table_func(section: str, func: int) -> int:
    # Proper dihedrals of types 1 and 9 live in the same table.
    return 1 if section == "dihedraltypes" and func == 9 else func


def _matches(pattern: tuple[str, ...], types: tuple[str, ...]) -> int:
    """Number of explicitly matched types, or -1 if the pattern does not match."""
    if all(p in (WILDCARD, t) for p, t in zip(pattern, types)):
        return sum(p != WILDCARD for p in pattern)
    return -1


class TypeIndex:
    """Lookup of the entries of the parameter tables of one topology."""

    def __init__(self, tables: dict[str, Iterable[TypeEntry]]) -> None:
        # (section, function) -> list of blocks of consecutive identical entries
        blocks: dict[tuple[str, int], list[list[TypeEntry]]] = defaultdict(list)
        slot: dict[tuple[str, int, tuple[str, ...]], int] = {}

        for section, entries in tables.items():
            for entry in entries:
                func = _table_func(section, entry.func)
                key = (section, func)
                same_block = (
                    section == "dihedraltypes"
                    and entry.func == 9
                    and blocks[key]
                    and blocks[key][-1][-1].func == 9
                    and blocks[key][-1][-1].types == entry.types
                )
                if same_block:
                    blocks[key][-1].append(entry)
                    continue
                slot_key = (*key, _canonical(entry.types))
                if slot_key in slot and entry.func != 9:
                    blocks[key][slot[slot_key]] = [entry]
                else:
                    slot[slot_key] = len(blocks[key])
                    blocks[key].append([entry])

        self._exact: dict[tuple[str, int, tuple[str, ...]], list[TypeEntry]] = {}
        self._wild: dict[tuple[str, int], list[list[TypeEntry]]] = defaultdict(list)
        for (section, func), block_list in blocks.items():
            for block in block_list:
                types = block[0].types
                if WILDCARD in types:
                    self._wild[section, func].append(block)
                else:
                    self._exact[section, func, _canonical(types)] = block

    def find(self, section: str, types: tuple[str, ...], func: int) -> list[TypeEntry]:
        """The table entries that apply to ``types``; empty if there is none."""
        func = _table_func(section, func)
        block = self._exact.get((section, func, _canonical(types)))
        if block is not None:
            return block

        best, best_score = None, -1
        for candidate in self._wild.get((section, func), ()):
            pattern = candidate[0].types
            score = max(_matches(pattern, types), _matches(pattern, types[::-1]))
            if score > best_score:
                best, best_score = candidate, score
        return best or []
