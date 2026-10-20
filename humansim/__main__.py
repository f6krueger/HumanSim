"""Command line: python -m humansim {show,coverage,check,compare,validate} ..."""

from __future__ import annotations

import argparse
import sys

from .checks import coupling_checks, growth_continuity, dataset_checks, modifier_checks, state_checks
from .context import Context, resolve
from .couplings import load_couplings
from .growth import load_growth
from .modifiers import load_modifiers
from .registry import Registry
from .solver import solve
from .report import format_compare, format_coverage, format_state, format_validation
from .sources import load_sources
from .units import parse_value
from .validation import load_validation


def _parse_edits(reg: Registry, items: list[str]) -> dict[str, float]:
    edits = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--set expects id=value, got {item!r}")
        pid, text = item.split("=", 1)
        pid = pid.strip()
        if pid not in reg:
            raise SystemExit(f"Unknown parameter {pid!r}. Run 'show --all' to list parameters.")
        edits[pid] = parse_value(text, reg[pid].unit)
    return edits


def _context(args, reg: Registry) -> Context:
    age: object = args.age
    try:
        age = float(args.age)
        age = int(age) if age.is_integer() else age
    except ValueError:
        pass
    return Context(sex=args.sex, age=age, growth_z=getattr(args, 'z', 0.0), edits=_parse_edits(reg, args.set),
                   keep_unadjusted=getattr(args, 'keep_unadjusted', False))


def _print_checks(results) -> int:
    worst = 0
    for r in results:
        print(f"[{r.status.upper():<7}] {r.id}: {r.message}")
        worst = max(worst, 1 if r.status == "fail" else 0)
    return worst


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="humansim", description="Reference human model")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_context_args(p):
        p.add_argument("--sex", default="male", choices=["male", "female"])
        p.add_argument("--age", default="adult", help="'adult' or years (below 20: CDC growth curves)")
        p.add_argument("--z", type=float, default=0.0,
                       help="pediatric ages only: growth-curve z-score (0 = median; 1 = about the 84th percentile)")
        p.add_argument("--set", action="append", metavar="ID=VALUE",
                       help="override a parameter, e.g. body_mass=90 or body_mass='198 lb'")
        p.add_argument("--solve", action="store_true",
                       help="run the steady-state solver so edits ripple through the couplings")
        p.add_argument("--keep-unadjusted", action="store_true",
                       help="keep other-sex values that have no sex rule (flagged no_sex_rule)")

    p_show = sub.add_parser("show", help="print resolved values")
    add_context_args(p_show)
    p_show.add_argument("--group")
    p_show.add_argument("--changed", action="store_true", help="only values that differ from baseline")
    p_show.add_argument("--all", action="store_true", help="include parameters with no value")
    p_show.add_argument("--trace", action="store_true", help="print the modifier/edit trace")

    p_cov = sub.add_parser("coverage", help="how many parameters have values and sources")
    p_cov.add_argument("--missing", action="store_true", help="list parameters with no value")

    p_chk = sub.add_parser("check", help="run consistency checks")
    add_context_args(p_chk)

    sub.add_parser("compare", help="male and female side by side")
    sub.add_parser("validate", help="sex rules vs stored pairs, and ANSUR II ratios vs an independent dataset")

    args = ap.parse_args(argv)
    reg = Registry.load()

    if args.cmd == "coverage":
        print(format_coverage(reg, list_missing=args.missing))
        return 0

    modifiers = load_modifiers()
    if args.cmd == "compare":
        print(format_compare(resolve(reg, Context(sex="male"), modifiers),
                             resolve(reg, Context(sex="female"), modifiers)))
        return 0
    if args.cmd == "validate":
        text, worst = format_validation(reg, modifiers, load_validation("ansur2_female"),
                                        load_validation("nchs_sex_ratios"))
        print(text)
        return 0 if worst < 1e-3 else 1
    try:
        ctx = _context(args, reg)
        state = (solve(reg, ctx, modifiers, load_couplings(reg)) if args.solve
                 else resolve(reg, ctx, modifiers))
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"error: {exc.args[0] if exc.args else exc}") from None

    if args.cmd == "show":
        print(format_state(state, group=args.group, only_changed=args.changed,
                           hide_missing=not args.all))
        for note in state.notes:
            print(f"\nNote: {note}")
        if state.solver:
            sv = state.solver
            print(f"\nSolver: {'converged' if sv['converged'] else 'DID NOT CONVERGE'} in {sv['iterations']} "
                  f"iteration(s); {len(sv['skipped'])} of {sv['couplings']} couplings skipped "
                  f"(missing data or pinned by an edit).")
        if args.trace:
            print("\nTrace:")
            for t in state.trace:
                b = "-" if t.before is None else f"{t.before:.6g}"
                a = "-" if t.after is None else f"{t.after:.6g}"
                print(f"  {t.stage:<8} {t.source:<28} {t.param}: {b} -> {a}  {t.note}")
        return 0

    sources = load_sources()
    couplings = load_couplings(reg)
    cited = ({m.source_ref for m in modifiers} | {c.source_ref for c in couplings}
             | {load_validation("nchs_sex_ratios")["source_ref"]}
             | {c['source_ref'] for c in load_growth()._curves.values()}) - {None}
    results = (state_checks(state) + dataset_checks(reg, sources, cited)
               + modifier_checks(reg, modifiers, sources)
               + coupling_checks(reg, couplings, sources)
               + growth_continuity(reg, state.context.sex))
    return _print_checks(results)


if __name__ == "__main__":
    sys.exit(main())
