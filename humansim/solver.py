"""Steady-state solver.

solve() resolves the context twice: once as given (with user edits) and once as the reference
(same sex and age, no edits). It then applies the couplings until the values stop moving.

* A parameter the user edited, or that a hormone, genotype or drug effect set, is pinned: couplings never overwrite it.
* Every other target is set to reference_value * g, where g is the coupling's ratio (see
  couplings.py). Compensation couplings have a deadband, so a change inside normal variation moves
  nothing; accounting couplings apply immediately.
* Couplings are applied in list order and repeated until the largest relative change is below
  `tol`. Feedback loops (a larger heart raises metabolic rate, which raises cardiac output, which
  asks for a larger heart) converge when their loop gain is below 1; otherwise the result says
  `converged: False` and the values must not be trusted.
* Missing inputs or a missing reference value skip the coupling and the skip is reported.
"""

from __future__ import annotations

import math
from typing import Iterable

from .context import Context, State, TraceEntry, finalize_derived, flag_ranges, resolve
from .couplings import FUNCTIONS, Coupling, SolveCtx
from .modifiers import Modifier
from .registry import Registry


def solve(registry: Registry, context: Context, modifiers: Iterable[Modifier],
          couplings: Iterable[Coupling], tol: float = 1e-9, max_iter: int = 200) -> State:
    modifiers, couplings = list(modifiers), list(couplings)
    # The reference is the same person without anything the user imposed: no edits, and no hormone
    # level, genotype or drug. Their effects therefore count as changes the couplings must respond to.
    ref = resolve(registry, Context(sex=context.sex, age=context.age, cycle_phase=context.cycle_phase,
                                    growth_z=context.growth_z, keep_unadjusted=context.keep_unadjusted), modifiers)
    state = resolve(registry, context, modifiers)
    for c in couplings:
        if c.target not in registry or registry[c.target].kind == "derived":
            raise KeyError(f"Coupling {c.id!r} targets unknown or derived parameter {c.target!r}")

    # Parameters the user edited, or that a hormone/genotype/drug effect set directly, are pinned.
    pinned = set(context.edits) | {t.param for t in state.trace
                                    if t.stage in ("hormone", "genotype", "drug") and t.after is not None}
    values = dict(state.values)
    # Derived values are recomputed from the current stored values inside each coupling that needs them.
    skipped: dict[str, str] = {}
    last: dict[str, tuple[float, float, str | None]] = {}  # coupling id -> (g, raw, ...)
    converged, iterations, diverged = False, 0, False
    for iterations in range(1, max_iter + 1):
        worst = 0.0
        for c in couplings:
            if c.target in pinned:
                skipped[c.id] = "target was edited by the user"
                continue
            base = ref.values.get(c.target)
            if base is None:
                skipped[c.id] = "no reference value for the target"
                continue
            ctx = SolveCtx(registry, context.sex, values, ref.values, pinned, c, state.ranges)
            try:
                out = FUNCTIONS[c.fn](ctx, c.args)
            except (OverflowError, ZeroDivisionError):
                diverged = True  # a runaway feedback loop: values left the representable range
                break
            if out is None:
                skipped[c.id] = "an input has no value"
                continue
            skipped.pop(c.id, None)
            g, raw = out
            new = base * g
            old = values[c.target]
            if old is not None:
                worst = max(worst, abs(new - old) / max(abs(old), 1e-12))
            values[c.target] = new
            last[c.id] = (g, raw)
        if diverged or not math.isfinite(worst):
            break
        if worst < tol:
            converged = True
            break

    flags = {k: set(v) for k, v in state.flags.items()}
    trace = list(state.trace)
    for c in couplings:
        if c.id not in last:
            continue
        g, raw = last[c.id]
        before, after = ref.values.get(c.target), values[c.target]
        if before is None or math.isclose(before, after, rel_tol=1e-9):
            continue
        flags[c.target].add("solved")
        if c.status != "referenced":
            flags[c.target].add(f"coupling_{c.status}")
        note = f"{c.kind} {c.role}; driving ratio {raw:.3g}, applied {g:.3g}"
        if c.validity and not (c.validity[0] <= raw <= c.validity[1]):
            flags[c.target].add("outside_coupling_validity")
            note += f"; OUTSIDE validity {c.validity[0]}-{c.validity[1]}: extrapolation"
        trace.append(TraceEntry(c.target, "solve", c.id, before, after, note))

    flag_ranges(registry, context.sex, values, flags, state.ranges, state.phase)
    missing_inputs = finalize_derived(values, flags)
    return State(context, registry, values, flags, trace, missing_inputs, ranges=state.ranges, notes=state.notes,
                 solver={"iterations": iterations, "converged": converged, "skipped": skipped,
                         "couplings": len(couplings)})
