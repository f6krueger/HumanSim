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
    # Keep values taken from the other sex that have no sex rule (flagged "no_sex_rule") instead of
    # leaving them empty. Off by default because they are not valid for this sex.
    keep_unadjusted: bool = False


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
    solver: dict | None = None  # filled by solver.solve(): iterations, converged, skipped couplings

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


def flag_ranges(registry: Registry, sex: str, values: dict, flags: dict) -> None:
    """Flag stored values outside their per-sex reference range or their viability limit."""
    for p in registry:
        v = values.get(p.id)
        if p.kind == "derived" or v is None:
            continue
        flags[p.id].discard("outside_reference_range")
        flags[p.id].discard("beyond_viability_limit")
        rng = p.ref_range.get(sex)  # ranges are per sex; no range for a sex means no flag
        if rng and not (rng[0] <= v <= rng[1]):
            flags[p.id].add("outside_reference_range")
        if p.viability_limit and not (p.viability_limit[0] <= v <= p.viability_limit[1]):
            flags[p.id].add("beyond_viability_limit")


def finalize_derived(values: dict, flags: dict) -> dict[str, list[str]]:
    """Compute derived values in place and set their flags. Returns the missing inputs per derived id."""
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
    return missing_inputs


def resolve(registry: Registry, context: Context, modifiers: Iterable[Modifier] = ()) -> State:
    if context.sex not in ("male", "female"):
        raise ValueError(f"sex must be 'male' or 'female', got {context.sex!r}")

    values: dict[str, float | None] = {}
    flags: dict[str, set[str]] = {}
    for p in registry:
        if p.kind == "derived":
            continue
        v = p.baseline(context.sex, context.age)
        f: set[str] = set()
        if v is None:
            # No stored value for this sex: start from the other sex's value and let a sex rule
            # convert it (forward for female, inverted for male). Stored values always win.
            other = "male" if context.sex == "female" else "female"
            ov = p.baseline(other, context.age)
            if ov is not None:
                v = ov
                f |= {f"derived_from_{other}", "no_sex_rule"}
        values[p.id] = v
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
    if context.sex == "male":
        # Sex rules are written female = f(male). A male value derived from a female-only value
        # uses the inverse of that rule (ratio -> divide, offset -> subtract; 'set' cannot be inverted).
        inverse = [m for m in modifiers if m.stage == "sex" and dict(m.when) == {"sex": "female"}]
        applicable = inverse + applicable
    derived_flag = f"derived_from_{'male' if context.sex == 'female' else 'female'}"
    for m in applicable:
        inverted = (m.stage == "sex" and context.sex == "male"
                    and dict(m.when) == {"sex": "female"})
        for pid, x in m.effects.items():
            if pid not in values:
                raise KeyError(f"Modifier {m.id!r} targets unknown or derived parameter {pid!r}")
            before = values[pid]
            if (m.stage == "sex" and (inverted or context.sex == "female")
                    and before is not None and derived_flag not in flags[pid]):
                # A stored value for this sex already is the answer; a rule would apply twice.
                if not inverted:
                    trace.append(TraceEntry(pid, m.stage, m.id, before, before,
                                            f"skipped: a stored {context.sex} value exists"))
                continue
            if before is None and (m.op != "set" or inverted):
                if not inverted:
                    flags[pid].add("modifier_skipped_no_baseline")
                    trace.append(TraceEntry(pid, m.stage, m.id, None, None, "skipped: no baseline value"))
                continue
            if inverted:
                if m.op == "set":
                    continue  # not invertible: the value stays flagged no_sex_rule
                after = before / x if m.op == "ratio" else before - x
            else:
                after = m.apply(before, x)
            values[pid] = after
            flags[pid].discard("missing")
            flags[pid].discard("no_data_for_context")
            if m.stage == "sex":
                flags[pid].discard("no_sex_rule")
            flags[pid].add("modified")
            if m.status != "referenced":
                flags[pid].add(f"modifier_{m.status}")
            trace.append(TraceEntry(pid, m.stage, m.id, before, after,
                                    ("inverse of female rule. " if inverted else "") + m.note))

    if not context.keep_unadjusted:
        for pid, f in flags.items():
            if "no_sex_rule" in f and values.get(pid) is not None:
                values[pid] = None  # a value from the other sex would be wrong for this body
                flags[pid] = {"missing", "no_sex_rule"}

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

    flag_ranges(registry, context.sex, values, flags)
    missing_inputs = finalize_derived(values, flags)

    return State(context, registry, values, flags, trace, missing_inputs)
