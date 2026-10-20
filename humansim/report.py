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


def format_validation(registry: Registry, modifiers, doc: dict, independent: dict | None = None) -> tuple[str, float]:
    """Two checks on the sex data. Returns (text, worst disagreement between a rule and its stored pair).

    1. Consistency: every sex rule must equal the ratio of the stored female and male values.
    2. Independent: female/male ratios from a second dataset (NCHS) against the ANSUR II ratios.
    """
    from .validation import scaling_analysis

    analysis = {r["id"]: r for r in scaling_analysis(doc, registry)}
    rules = {pid: m for m in modifiers if m.stage == "sex" and m.op == "ratio" for pid in m.effects}
    lines = ["1. Consistency: sex rule ratio vs the stored female/male pair.",
             "   Both sexes are stored; the ratios are computed from those pairs, so they agree by",
             "   construction. This catches stale rule files, not scientific error.", "",
             f"{'parameter':<26}{'male':>10}{'female':>10}{'F/M':>8}{'rule':>8}  rule type"]
    worst = 0.0
    for p in registry:
        m, f = p.baseline("male", "adult"), p.baseline("female", "adult")
        if p.kind == "derived" or m is None or f is None:
            continue
        rule = rules.get(p.id)
        rr = rule.effects[p.id] if rule else None
        if rr is None:
            lines.append(f"{p.id:<26}{m:>10.2f}{f:>10.2f}{f / m:>8.3f}{'-':>8}  NO RULE")
            worst = float("inf")
            continue
        worst = max(worst, abs(rr - f / m))
        lines.append(f"{p.id:<26}{m:>10.2f}{f:>10.2f}{f / m:>8.3f}{rr:>8.3f}  {rule.rule_type}")
    lines += ["", f"Largest rule/pair disagreement: {worst:.5f}"]

    if independent:
        lines += ["", f"2. Independent check: {independent.get('population', 'second dataset')}.",
                  "   Not used to build any rule. A different population (general public vs US Army),",
                  "   so differences are expected; large ones would flag a problem.", "",
                  f"{'parameter':<26}{'model F/M':>10}{'check F/M':>10}{'diff':>8}"]
        for pid, rec in independent["parameters"].items():
            m, f = registry[pid].baseline("male", "adult"), registry[pid].baseline("female", "adult")
            if m is None or f is None:
                continue
            model, check = f / m, rec["female"] / rec["male"]
            lines.append(f"{pid:<26}{model:>10.3f}{check:>10.3f}{model / check - 1:>+8.1%}")
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
