"""Modifiers: the single mechanism for sex, hormone, genotype and drug effects.

A modifier changes one or more parameters by ratio, offset or replacement, and applies
only when its `when` conditions match the Context. Stages run in a fixed order
(sex, hormone, genotype, drug) and user edits are applied after all of them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import yaml

from .registry import DATA_DIR, RegistryError

STAGES = ("sex", "hormone", "genotype", "drug")
OPS = ("ratio", "offset", "set")


@dataclass(frozen=True)
class Modifier:
    id: str
    stage: str
    op: str
    effects: Mapping[str, float]
    when: Mapping[str, object] = field(default_factory=dict)
    status: str = "unsourced"
    source_ref: str | None = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.stage not in STAGES:
            raise RegistryError(f"Modifier {self.id!r}: stage must be one of {STAGES}")
        if self.op not in OPS:
            raise RegistryError(f"Modifier {self.id!r}: op must be one of {OPS}")
        if not self.effects:
            raise RegistryError(f"Modifier {self.id!r}: no effects")

    def applies(self, ctx) -> bool:
        """`when` keys: sex, age (equality); hormone:<name>, drug:<name> (present in context);
        genotype:<gene> (equals the given allele string)."""
        for key, want in self.when.items():
            if key == "sex":
                ok = ctx.sex == want
            elif key == "age":
                ok = ctx.age == want
            elif key.startswith("hormone:"):
                ok = key.split(":", 1)[1] in ctx.hormone_state
            elif key.startswith("drug:"):
                ok = key.split(":", 1)[1] in ctx.drugs
            elif key.startswith("genotype:"):
                ok = ctx.genotype.get(key.split(":", 1)[1]) == want
            else:
                raise RegistryError(f"Modifier {self.id!r}: unknown condition {key!r}")
            if not ok:
                return False
        return True

    def apply(self, before: float | None, x: float) -> float:
        if self.op == "set":
            return float(x)
        assert before is not None
        return before * x if self.op == "ratio" else before + x


def load_modifiers(data_dir: str | Path | None = None) -> list[Modifier]:
    root = Path(data_dir) if data_dir else DATA_DIR
    out: list[Modifier] = []
    for path in sorted((root / "modifiers").glob("*.yaml")):
        doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for raw in doc.get("modifiers") or []:
            out.append(Modifier(**raw))
    return out
