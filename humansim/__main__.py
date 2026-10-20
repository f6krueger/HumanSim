"""Command line: python -m humansim {show,coverage,check,compare,validate} ..."""

from __future__ import annotations

import argparse
import sys

from .checks import dataset_checks, modifier_checks, state_checks
from .context import Context, resolve
from .modifiers import load_modifiers
from .registry import Registry
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
    return Context(sex=args.sex, age=age, edits=_parse_edits(reg, args.set),
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
        p.add_argument("--age", default="adult")
        p.add_argument("--set", action="append", metavar="ID=VALUE",
                       help="override a parameter, e.g. body_mass=90 or body_mass='198 lb'")
        p.add_argument("--keep-unadjusted", action="store_true",
                       help="female only: keep male values that have no sex rule (flagged no_sex_rule)")

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

    sub.add_parser("compare", help="male and derived female side by side")
    sub.add_parser("validate", help="derived female vs the ANSUR II female means")

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
        text, worst = format_validation(reg, modifiers, load_validation("ansur2_female"))
        print(text)
        return 0 if worst < 1e-3 else 1
    try:
        state = resolve(reg, _context(args, reg), modifiers)
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"error: {exc.args[0] if exc.args else exc}") from None

    if args.cmd == "show":
        print(format_state(state, group=args.group, only_changed=args.changed,
                           hide_missing=not args.all))
        if args.trace:
            print("\nTrace:")
            for t in state.trace:
                b = "-" if t.before is None else f"{t.before:.6g}"
                a = "-" if t.after is None else f"{t.after:.6g}"
                print(f"  {t.stage:<8} {t.source:<28} {t.param}: {b} -> {a}  {t.note}")
        return 0

    sources = load_sources()
    results = (state_checks(state) + dataset_checks(reg, sources)
               + modifier_checks(reg, modifiers, sources))
    return _print_checks(results)


if __name__ == "__main__":
    sys.exit(main())
