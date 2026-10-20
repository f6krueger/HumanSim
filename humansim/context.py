"""Context and resolution.

resolve() turns a Context into a State:
  1. baseline values for (sex, age) from the registry
  2. modifiers, in stage order (sex, hormone, genotype, drug)
  3. user edits
  4. range checks on stored parameters
  5. derived values
Every change is recorded in the trace, and every value carries flags.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

from . import derived as _derived
from .modifiers import STAGES, Modifier
from .registry import Registry

# Flags a derived value inherits from its inputs.
_INHERITED = {"unsourced", "placeholder", "edited", "modified"}


@dataclass
class Context:
    sex: str = "male"
    age: object = "adult"
    hormone_state: dict[str, float] = field(default_factory=dict)
    genotype: dict[str, str] = field(default_factory=dict)
    drugs: dict[str, float] = field(default_factory=dict)
    edits: dict[str, float] = field(default_factory=dict)  # parameter id -> value in registry unit


@dataclass(frozen=True)
class TraceEntry:
    param: str
    stage: str  # sex | hormone | genotype | drug | edit
    source: str  # modifier id, or "user"
    before: float | None
    after: float | None
    note: str = ""


@dataclass
class State:
    context: Context
    registry: Registry
    values: dict[str, float | None]
    flags: dict[str, set[str]]
    trace: list[TraceEntry]
    missing_inputs: dict[str, list[str]]

    def get(self, pid: str) -> float | None:
        return self.values[pid]

    def changed_from(self, other: "State", rel_tol: float = 1e-9) -> dict[str, tuple[float | None, float | None]]:
        out = {}
        for pid, v in self.values.items():
            w = other.values.get(pid)
            if v is None and w is None:
                continue
            if v is None or w is None or not math.isclose(v, w, rel_tol=rel_tol, abs_tol=1e-12):
                out[pid] = (w, v)
        return out


def resolve(registry: Registry, context: Context, modifiers: Iterable[Modifier] = ()) -> State:
    if context.sex not in ("male", "female"):
        raise ValueError(f"sex must be 'male' or 'female', got {context.sex!r}")

    values: dict[str, float | None] = {}
    flags: dict[str, set[str]] = {}
    for p in registry:
        if p.kind == "derived":
            continue
        v = p.baseline(context.sex, context.age)
        values[p.id] = v
        f: set[str] = set()
        if v is None:
            f.add("missing")
            if p.values:
                f.add("no_data_for_context")
        elif p.status in ("unsourced", "placeholder"):
            f.add(p.status)
        flags[p.id] = f

    trace: list[TraceEntry] = []

    applicable = [m for m in modifiers if m.applies(context)]
    applicable.sort(key=lambda m: STAGES.index(m.stage))  # stable within a stage
    for m in applicable:
        for pid, x in m.effects.items():
            if pid not in values:
                raise KeyError(f"Modifier {m.id!r} targets unknown or derived parameter {pid!r}")
            before = values[pid]
            if before is None and m.op != "set":
                flags[pid].add("modifier_skipped_no_baseline")
                trace.append(TraceEntry(pid, m.stage, m.id, None, None, "skipped: no baseline value"))
                continue
            after = m.apply(before, x)
            values[pid] = after
            flags[pid].discard("missing")
            flags[pid].discard("no_data_for_context")
            flags[pid].add("modified")
            if m.status != "referenced":
                flags[pid].add(f"modifier_{m.status}")
            trace.append(TraceEntry(pid, m.stage, m.id, before, after, m.note))

    for pid, x in context.edits.items():
        if pid not in registry:
            raise KeyError(f"Unknown parameter {pid!r}")
        if registry[pid].kind == "derived":
            raise ValueError(f"{pid!r} is derived and cannot be edited; edit one of its inputs")
        x = float(x)
        if not math.isfinite(x):
            raise ValueError(f"Edit for {pid!r} must be a finite number")
        before = values[pid]
        values[pid] = x
        flags[pid] = {"edited"}
        trace.append(TraceEntry(pid, "edit", "user", before, x))

    for p in registry:
        v = values.get(p.id)
        if p.kind == "derived" or v is None:
            continue
        if p.ref_range and not (p.ref_range[0] <= v <= p.ref_range[1]):
            flags[p.id].add("outside_reference_range")
        if p.viability_limit and not (p.viability_limit[0] <= v <= p.viability_limit[1]):
            flags[p.id].add("beyond_viability_limit")

    missing_inputs: dict[str, list[str]] = {}
    dvals = _derived.compute(values)
    values.update(dvals)
    for _ in range(3):  # a few passes cover derived-of-derived chains
        for spec in _derived.SPECS:
            gone = [i for i in spec.inputs if values.get(i) is None]
            f = {"derived"}
            if dvals[spec.id] is None:
                f.add("missing")
                missing_inputs[spec.id] = gone
            for i in spec.inputs:
                f |= flags.get(i, set()) & _INHERITED
            flags[spec.id] = f

    return State(context, registry, values, flags, trace, missing_inputs)
