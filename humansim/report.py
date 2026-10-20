"""Text reports: data coverage and resolved state."""

from __future__ import annotations

from collections import Counter

from .context import State
from .registry import Registry


def coverage_counts(registry: Registry) -> dict[str, Counter]:
    out: dict[str, Counter] = {}
    for p in registry:
        out.setdefault(p.group, Counter())[p.status] += 1
    return out


def format_coverage(registry: Registry, list_missing: bool = False) -> str:
    cols = ["referenced", "unsourced", "placeholder", "missing", "derived"]
    lines = [f"{'group':<16}" + "".join(f"{c:>12}" for c in cols)]
    total: Counter = Counter()
    for group, counts in coverage_counts(registry).items():
        lines.append(f"{group:<16}" + "".join(f"{counts[c]:>12}" for c in cols))
        total.update(counts)
    lines.append(f"{'TOTAL':<16}" + "".join(f"{total[c]:>12}" for c in cols))
    stored = sum(v for k, v in total.items() if k != "derived")
    have = total["referenced"] + total["unsourced"] + total["placeholder"]
    lines.append("")
    lines.append(f"{have} of {stored} stored parameters have a value; "
                 f"{total['referenced']} of them have a source.")
    if list_missing:
        lines.append("")
        lines.append("Missing values:")
        for p in registry:
            if p.status == "missing":
                lines.append(f"  {p.id:<28} {p.unit:<14} {p.group}")
    return "\n".join(lines)


def _fmt(v: float | None) -> str:
    return "-" if v is None else f"{v:.4g}"


def format_state(state: State, group: str | None = None, only_changed: bool = False,
                 hide_missing: bool = True) -> str:
    reg = state.registry
    base = None
    if only_changed:
        from .context import Context, resolve
        base = resolve(reg, Context(sex=state.context.sex, age=state.context.age))
        changed = state.changed_from(base)
    rows = []
    for p in reg:
        if group and p.group != group:
            continue
        v = state.get(p.id)
        if hide_missing and v is None and not only_changed:
            continue
        if only_changed and p.id not in changed:
            continue
        flags = ",".join(sorted(state.flags.get(p.id, ())))
        shown = _fmt(v)
        if only_changed:
            shown = f"{_fmt(changed[p.id][0])} -> {_fmt(v)}"
        rows.append((p.id, shown, p.unit, flags))
    if not rows:
        return "(nothing to show)"
    w0 = max(len(r[0]) for r in rows)
    w1 = max(len(r[1]) for r in rows)
    w2 = max(len(r[2]) for r in rows)
    return "\n".join(f"{a:<{w0}}  {b:>{w1}}  {c:<{w2}}  {d}" for a, b, c, d in rows)
