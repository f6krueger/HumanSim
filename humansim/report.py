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


def format_compare(male: State, female: State) -> str:
    rows = []
    for p in male.registry:
        m, f = male.get(p.id), female.get(p.id)
        if m is None and f is None:
            continue
        ratio = "-" if (m is None or f is None or m == 0) else f"{f / m:.3f}"
        notes = ",".join(sorted(female.flags.get(p.id, set()) & {"no_sex_rule", "unsourced", "missing"}))
        rows.append((p.id, _fmt(m), _fmt(f), ratio, p.unit, notes))
    w = [max(len(r[i]) for r in rows) for i in range(5)]
    head = f"{'parameter':<{w[0]}}  {'male':>{w[1]}}  {'female':>{w[2]}}  {'F/M':>{w[3]}}  {'unit':<{w[4]}}"
    return head + "\n" + "\n".join(
        f"{a:<{w[0]}}  {b:>{w[1]}}  {c:>{w[2]}}  {d:>{w[3]}}  {e:<{w[4]}}  {n}" for a, b, c, d, e, n in rows)


def format_validation(registry: Registry, modifiers, doc: dict) -> tuple[str, float]:
    """Compare the derived female with the ANSUR II female means. Returns (text, worst % difference)."""
    from .context import Context, resolve
    from .validation import scaling_analysis

    female = resolve(registry, Context(sex="female"), modifiers)
    analysis = {r["id"]: r for r in scaling_analysis(doc, registry)}
    lines = [
        "Calibration check: derived female vs ANSUR II female means.",
        "The rules were fitted to these means, so agreement is expected. This is a consistency",
        "check, NOT independent validation.",
        "",
        f"{'parameter':<26}{'derived':>10}{'ANSUR II':>10}{'diff':>8}  rule type",
    ]
    worst = 0.0
    for pid, rec in doc["parameters"].items():
        d = female.get(pid)
        if d is None:
            lines.append(f"{pid:<26}{'-':>10}{rec['mean']:>10.2f}{'':>8}  NO DERIVED VALUE")
            worst = float("inf")
            continue
        diff = d / rec["mean"] - 1.0
        worst = max(worst, abs(diff))
        lines.append(f"{pid:<26}{d:>10.2f}{rec['mean']:>10.2f}{diff:>+8.2%}  {analysis[pid]['type']}")
    lines.append("")
    lines.append(f"Largest difference: {worst:.3%}")
    return "\n".join(lines), worst


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
