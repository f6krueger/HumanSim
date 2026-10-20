"""Command line: python -m humansim {show,coverage,check,compare,validate} ..."""

from __future__ import annotations

import argparse
import sys

from .checks import drug_checks, coupling_checks, growth_continuity, dataset_checks, modifier_checks, state_checks
from .context import Context, resolve
from .couplings import load_couplings
from .pharmacology import Dose, load_drugs, simulate as pk_simulate
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
    hormones = {}
    for item in getattr(args, "hormone", None) or []:
        pid, _, text = item.partition("=")
        if pid not in reg or reg[pid].group != "hormones" or not text:
            raise SystemExit(f"--hormone needs ID=VALUE for a hormone parameter, got {item!r}")
        hormones[pid] = parse_value(text, reg[pid].unit)
    kwargs = dict(sex=args.sex, age=age, growth_z=getattr(args, "z", 0.0), edits=_parse_edits(reg, args.set),
                  keep_unadjusted=getattr(args, "keep_unadjusted", False), hormone_state=hormones,
                  cycle_phase=getattr(args, "cycle_phase", None))
    day = getattr(args, "cycle_day", None)
    if day is not None:
        if args.sex != "female":
            raise SystemExit("--cycle-day applies to --sex female")
        from .cycle import cycle_context
        ctx = cycle_context(reg, day, age=age, edits=kwargs["edits"], keep_unadjusted=kwargs["keep_unadjusted"])
        ctx.hormone_state.update(hormones)
        return ctx
    return Context(**kwargs)


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
        p.add_argument("--hormone", action="append", metavar="ID=VALUE",
                       help="set a hormone level, e.g. testosterone_total=1345 (ng/dL); drives hormone effects")
        p.add_argument("--cycle-phase", choices=["follicular", "midcycle", "luteal", "postmenopausal"],
                       help="adult women: use this phase's reference hormone values")
        p.add_argument("--cycle-day", type=float, help="adult women: cycle day 1-28 (interpolated profile)")
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

    p_pk = sub.add_parser("pk", help="drug concentration over time after intravenous doses")
    add_context_args(p_pk)
    p_pk.add_argument("--drug", required=True)
    p_pk.add_argument("--dose", type=float, required=True, help="mg")
    p_pk.add_argument("--infusion", type=float, default=0.0, help="infusion time in hours (0 = bolus)")
    p_pk.add_argument("--hours", type=float, default=24.0)

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

    if args.cmd == "pk":
        drugs = load_drugs()
        if args.drug not in drugs:
            raise SystemExit(f"error: unknown drug {args.drug!r}; known: {', '.join(drugs)}")
        try:
            r = pk_simulate(reg, state, drugs[args.drug], [Dose(args.dose, 0.0, args.infusion)], args.hours,
                            modifiers)
        except ValueError as exc:
            raise SystemExit(f"error: {exc}") from None
        print(f"{r.drug.name}, {args.dose:g} mg IV, {state.context.sex}")
        print(f"clearance {r.clearance_l_h:.1f} L/h, volume {r.volume_l:.1f} L, half-life {r.half_life_h:.2f} h")
        print(f"Cmax {r.cmax_mg_l:.2f} mg/L, AUC over {args.hours:g} h {r.auc_mg_h_l():.1f} mg*h/L")
        print(f"flags: {', '.join(sorted(r.flags))}")
        print(f"note: {r.notes[0]}")
        return 0

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
             | {load_validation("nchs_sex_ratios")["source_ref"], "mifflin1990"}
             | {c['source_ref'] for c in load_growth()._curves.values()}
             | {d.source_ref for d in load_drugs().values()}) - {None}
    results = (state_checks(state) + dataset_checks(reg, sources, cited)
               + modifier_checks(reg, modifiers, sources)
               + coupling_checks(reg, couplings, sources)
               + growth_continuity(reg, state.context.sex)
               + drug_checks(reg, load_drugs(), sources, modifiers))
    return _print_checks(results)


if __name__ == "__main__":
    sys.exit(main())
