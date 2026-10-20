"""Command line: python -m humansim {show,coverage,check,compare,validate} ..."""

from __future__ import annotations

import argparse
import sys

from .checks import drug_checks, coupling_checks, growth_continuity, dataset_checks, modifier_checks, state_checks
from .context import Context, resolve
from .couplings import load_couplings
from .genetics import load_genes, validate_genotype
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


def _parse_gene_edits(args) -> tuple:
    """--edit LOCUS:STATE (catalogue) and --hypothesis GENE:PARAM=OP:VALUE (your own idea)."""
    from .editing import Effect, GeneEdit, load_catalogue

    common = dict(mode=getattr(args, "edit_mode", "somatic"), fraction=getattr(args, "edit_fraction", 1.0),
                  tissue=getattr(args, "edit_tissue", None))
    out = []
    try:
        for item in getattr(args, "edit", None) or []:
            locus, _, state = item.partition(":")
            cat = load_catalogue()
            if locus not in cat or not state:
                raise SystemExit(f"--edit needs LOCUS:STATE from the catalogue ({', '.join(cat)}), got {item!r}")
            out.append(GeneEdit(gene=cat[locus].gene, locus=locus, state=state, **common))
        for item in getattr(args, "hypothesis", None) or []:
            gene, _, rest = item.partition(":")
            param, _, spec = rest.partition("=")
            op, _, val = spec.partition(":")
            if not (gene and param and op and val):
                raise SystemExit(f"--hypothesis needs GENE:PARAM=OP:VALUE (e.g. MSTN:skeletal_muscle_mass=ratio:1.2), got {item!r}")
            eff = Effect(param=param, op=op, value=float(val), tier="hypothetical",
                         rationale="supplied by the user as a hypothesis; not data")
            out.append(GeneEdit(gene=gene, hypothesis=(eff,), **common))
    except (ValueError, KeyError) as exc:
        raise SystemExit(f"error: {exc.args[0] if exc.args else exc}") from None
    return tuple(out)


def _catalogue_check(reg: Registry, sources: dict):
    from .checks import CheckResult
    from .editing import catalogue_checks, load_catalogue

    cat = load_catalogue()
    problems = catalogue_checks(reg, cat, sources)
    n_meas = sum(e.tier == "measured" for loc in cat.values() for _, effs in loc.states.values() for e in effs)
    n_hyp = sum(e.tier == "hypothetical" for loc in cat.values() for _, effs in loc.states.values() for e in effs)
    return [CheckResult("gene_catalogue", "fail" if problems else "pass",
                        "; ".join(problems) if problems else
                        f"{len(cat)} loci: {n_meas} measured, {n_hyp} hypothetical effects; all parameters and sources exist.")]


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
    genotype = {}
    for item in getattr(args, "genotype", None) or []:
        gid, _, dip = item.partition("=")
        if not dip:
            raise SystemExit(f"--genotype needs GENE=DIPLOTYPE like CYP2C9='*1/*3', got {item!r}")
        genotype[gid] = dip
    try:
        validate_genotype(load_genes(), genotype)
    except (KeyError, ValueError) as exc:
        raise SystemExit(f"error: {exc.args[0] if exc.args else exc}") from None
    gene_edits = _parse_gene_edits(args)
    kwargs = dict(sex=args.sex, genotype=genotype, gene_edits=gene_edits, age=age, growth_z=getattr(args, "z", 0.0), edits=_parse_edits(reg, args.set),
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
        p.add_argument("--genotype", action="append", metavar="GENE=DIPLOTYPE",
                       help="e.g. CYP2C9='*1/*3'; changes drug clearance where the drug has data")
        p.add_argument("--edit", action="append", metavar="LOCUS:STATE",
                       help="apply a catalogue gene edit, e.g. MSTN:loss_of_function (see 'catalogue')")
        p.add_argument("--hypothesis", action="append", metavar="GENE:PARAM=OP:VALUE",
                       help="your own hypothetical edit, e.g. MSTN:skeletal_muscle_mass=ratio:1.2 (always flagged hypothetical)")
        p.add_argument("--edit-mode", default="somatic", choices=["somatic", "germline"])
        p.add_argument("--edit-fraction", type=float, default=1.0, help="somatic: share of cells edited (0-1]")
        p.add_argument("--edit-tissue", help="somatic: tissue edited (recorded only)")
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

    p_genes = sub.add_parser("genes", help="pharmacogenes, activity scores and population frequencies")
    p_genes.add_argument("--population", default="european")

    p_serve = sub.add_parser("serve", help="start the local web UI")
    p_serve.add_argument("--port", type=int, default=8765)
    p_serve.add_argument("--open", action="store_true", help="open the page in your browser")

    sub.add_parser("catalogue", help="gene-edit catalogue with the tier of every effect")

    sub.add_parser("compare", help="male and female side by side")
    sub.add_parser("validate", help="sex rules vs stored pairs, and ANSUR II ratios vs an independent dataset")

    args = ap.parse_args(argv)
    if args.cmd == "serve":
        from .server import serve
        serve(args.port, args.open)
        return 0
    reg = Registry.load()

    if args.cmd == "coverage":
        print(format_coverage(reg, list_missing=args.missing))
        return 0

    modifiers = load_modifiers()
    if args.cmd == "catalogue":
        from .editing import load_catalogue, load_gene_index
        print("Every effect is labelled. HYPOTHETICAL means a stated guess, not data.\n")
        for loc in load_catalogue().values():
            print(f"{loc.id}: {loc.name}")
            for sname, (desc, effects) in loc.states.items():
                print(f"  state {sname}: {desc}")
                for e in effects:
                    tag = "MEASURED" if e.tier == "measured" else "HYPOTHETICAL"
                    print(f"    [{tag}] {e.param} {e.op} {e.value:g}")
        n = len(load_gene_index())
        print(f"\nGene index: {n} symbols" + ("" if n else " (none imported; see tools/import_hgnc.py). Any symbol is accepted for --hypothesis."))
        return 0
    if args.cmd == "genes":
        for g in load_genes().values():
            print(f"{g.id}: {g.name} (activity scores from {g.activity_source_ref})")
            try:
                dips = g.diplotype_frequencies(args.population)
            except KeyError as exc:
                print(f"  {exc.args[0]}")
                continue
            for dip, f in sorted(dips.items(), key=lambda kv: -kv[1]):
                print(f"  {dip:<8} score {g.activity_score(dip):<4g} {g.phenotype(dip):<13} {f:6.2%}")
            print("  by phenotype: " + ", ".join(f"{k} {v:.1%}" for k, v in g.phenotype_frequencies(args.population).items()))
        return 0
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
        for note in r.notes:
            print(f"note: {note}")
        return 0

    if args.cmd == "show":
        print(format_state(state, group=args.group, only_changed=args.changed,
                           hide_missing=not args.all))
        for note in state.notes:
            print(f"\nNote: {note}")
        if any("hypothetical" in f or "downstream_of_hypothetical" in f for f in state.flags.values()):
            print("\n*** HYPOTHETICAL: this state includes gene-edit effects that are NOT data (flagged hypothetical / "
                  "downstream_of_hypothetical). It is an exploration, not a prediction. ***")
        if any("germline" in f for f in state.flags.values()):
            print("Germline editing is heritable and prohibited or restricted for human embryos in many countries.")
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
             | {d.source_ref for d in load_drugs().values()}
             | {e['source_ref'] for d in load_drugs().values() for e in d.genotype_effects}
             | {x for g in load_genes().values() for x in [g.activity_source_ref] + [p['source_ref'] for p in g.populations.values()]}) - {None}
    results = (state_checks(state) + dataset_checks(reg, sources, cited)
               + modifier_checks(reg, modifiers, sources)
               + coupling_checks(reg, couplings, sources)
               + growth_continuity(reg, state.context.sex)
               + drug_checks(reg, load_drugs(), sources, modifiers)
               + _catalogue_check(reg, sources))
    return _print_checks(results)


if __name__ == "__main__":
    sys.exit(main())
