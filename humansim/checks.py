"""Consistency checks.

state_checks() inspects one resolved State. dataset_checks() inspects the registry
against the sources sidecar. A check with too little data reports "skipped", never "pass".
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .context import State
from .registry import Registry

# ICRP 89 Table 2.8: the rows that make up the whole-body sum (rows the report marks as duplicating
# others are left out). The report says they represent 96% of total body mass.
WHOLE_BODY_ROWS = (
    "separable_adipose_mass", "adrenals_mass", "tongue_mass", "salivary_glands_mass", "oesophagus_wall_mass",
    "stomach_wall_mass", "stomach_contents_mass", "small_intestine_wall_mass", "small_intestine_contents_mass",
    "right_colon_wall_mass", "right_colon_contents_mass", "left_colon_wall_mass", "left_colon_contents_mass",
    "rectosigmoid_wall_mass", "rectosigmoid_contents_mass", "liver_mass", "gallbladder_wall_mass",
    "gallbladder_contents_mass", "pancreas_mass", "brain_mass", "breasts_mass", "heart_mass", "blood_mass",
    "eyes_mass", "skin_mass", "skeletal_muscle_mass", "pituitary_mass", "larynx_mass", "trachea_mass",
    "lungs_mass", "cortical_bone_mass", "trabecular_bone_mass", "active_marrow_mass", "inactive_marrow_mass",
    "cartilage_mass", "teeth_mass", "skeleton_miscellaneous_mass", "spleen_mass", "thymus_mass", "thyroid_mass",
    "tonsils_mass", "kidneys_mass", "ureters_mass", "urinary_bladder_mass", "urethra_mass", "testes_mass",
    "epididymes_mass", "prostate_mass", "ovaries_mass", "fallopian_tubes_mass", "uterus_mass",
)
SEX_ORGANS = {"testes_mass", "epididymes_mass", "prostate_mass", "ovaries_mass", "fallopian_tubes_mass", "uterus_mass",
              "breasts_mass"}
ICRP_COVERAGE = 0.96  # footnote c of Table 2.8
COVERAGE_TOLERANCE = 0.03


@dataclass
class CheckResult:
    id: str
    status: str  # pass | warn | fail | skipped
    message: str
    details: dict = field(default_factory=dict)


def organ_reference_mass(state: State) -> tuple[float | None, str]:
    """The mass the organ data should add up to. The ICRP organ masses describe the ICRP reference person
    (73 kg man, 60 kg woman), not the ANSUR II body, so an unedited state is compared with the reference
    mass. Once something is edited or solved, body_mass (which includes the change) is used."""
    ref = state.get("reference_body_mass") if "reference_body_mass" in state.values else None
    if ref is not None and not state.context.edits and state.solver is None:
        return ref, "reference_body_mass"
    return state.get("body_mass"), "body_mass"


def tissue_mass_budget(state: State) -> CheckResult:
    mass, basis = organ_reference_mass(state)
    if mass is None:
        return CheckResult("tissue_mass_budget", "skipped", "no body mass for this context")
    sex_ok = [i for i in WHOLE_BODY_ROWS if i in state.values]
    missing, total = [], 0.0
    for i in sex_ok:
        v = state.get(i)
        if v is None:
            if i in SEX_ORGANS:
                continue  # not every sex or age has these organs
            missing.append(i)
        else:
            total += v
    details = {"listed_total_kg": round(total, 3), "basis": basis, "basis_mass_kg": round(mass, 3),
               "unassigned_kg": round(mass - total, 3), "missing": missing}
    if total > mass:
        return CheckResult("tissue_mass_budget", "fail",
                           f"Listed tissues sum to {total:.2f} kg, more than {basis} {mass:.2f} kg.", details)
    if missing:
        return CheckResult("tissue_mass_budget", "warn",
                           f"Partial budget: {len(missing)} tissue(s) have no value ({', '.join(missing[:4])}"
                           f"{'...' if len(missing) > 4 else ''}); {mass - total:.2f} kg of {basis} unassigned.", details)
    frac = total / mass
    if abs(frac - ICRP_COVERAGE) > COVERAGE_TOLERANCE:
        return CheckResult("tissue_mass_budget", "warn",
                           f"Listed tissues are {frac:.1%} of {basis} ({mass:.1f} kg). The ICRP 89 reference "
                           f"values cover 96%; the rest is connective and lymphatic tissue the report does not list.",
                           details)
    return CheckResult("tissue_mass_budget", "pass",
                       f"Listed tissues are {frac:.1%} of {basis} ({mass:.1f} kg); ICRP 89 reports 96%, with the rest "
                       f"separable connective and lymphatic tissue.", details)


def segment_stature_sum(state: State, warn: float = 0.05, fail: float = 0.10) -> CheckResult:
    parts = ("head_height", "neck_length", "torso_length", "crotch_height")
    vals = [state.get(i) for i in parts]
    height = state.get("height")
    if height is None or any(v is None for v in vals):
        return CheckResult("segment_stature_sum", "skipped", "segment values missing")
    total = sum(vals)
    err = abs(total - height) / height
    details = {"sum_cm": round(total, 2), "height_cm": height, "relative_error": round(err, 4)}
    if err > fail:
        return CheckResult("segment_stature_sum", "fail",
                           f"Head + neck + torso + crotch height = {total:.1f} cm vs stature {height:.1f} cm.",
                           details)
    if err > warn:
        return CheckResult("segment_stature_sum", "warn",
                           f"Segments differ from stature by {err:.1%}.", details)
    return CheckResult("segment_stature_sum", "pass",
                       f"Segments sum to within {err:.1%} of stature.", details)


def range_flags(state: State) -> CheckResult:
    outside = sorted(i for i, f in state.flags.items() if "outside_reference_range" in f)
    beyond = sorted(i for i, f in state.flags.items() if "beyond_viability_limit" in f)
    if beyond:
        return CheckResult("range_flags", "fail", f"Beyond viability limit: {', '.join(beyond)}",
                           {"beyond_viability": beyond, "outside_reference": outside})
    if outside:
        return CheckResult("range_flags", "warn", f"Outside reference range: {', '.join(outside)}",
                           {"outside_reference": outside})
    return CheckResult("range_flags", "pass", "All stored values with a known range are inside it.")


def sex_rule_coverage(state: State) -> CheckResult:
    """Parameters whose value exists only for the other sex and have no sex rule to convert it."""
    other = "male" if state.context.sex == "female" else "female"
    no_rule = sorted(i for i, f in state.flags.items() if "no_sex_rule" in f and i not in SEX_ORGANS)
    if no_rule:
        shown = ", ".join(no_rule[:8]) + (f", and {len(no_rule) - 8} more" if len(no_rule) > 8 else "")
        return CheckResult(
            "sex_rule_coverage", "warn",
            f"{len(no_rule)} parameter(s) have a {other} value but no {state.context.sex} value and no "
            f"sex rule, so they are empty (or flagged, with keep_unadjusted): {shown}",
            {"no_sex_rule": no_rule})
    return CheckResult("sex_rule_coverage", "pass",
                       f"Every parameter with a {other} value also has a {state.context.sex} value or a sex rule.")


RMR_RATES = {"liver_mass": 200.0, "brain_mass": 240.0, "heart_mass": 440.0, "kidneys_mass": 440.0,
             "skeletal_muscle_mass": 13.0, "adipose_tissue_mass": 4.5}
RMR_RESIDUAL_RATE = 12.0  # kcal/kg/day (Wang 2010, secondary source)


def rmr_organ_crosscheck(state: State, warn: float = 0.15) -> CheckResult:
    """Independent cross-check: the organ-sum metabolic rate (Wang 2010 tissue rates) against the ICRP
    basal metabolic rate. Two different sources, so agreement is a real check."""
    if state.ranges is not None:
        return CheckResult("rmr_organ_crosscheck", "skipped",
                           "the tissue rates (Wang 2010) are adult values; children's organs burn more per kg")
    mass, basis = organ_reference_mass(state)
    need = list(RMR_RATES) + ["resting_metabolic_rate"]
    if mass is None or any(state.get(i) is None for i in need):
        return CheckResult("rmr_organ_crosscheck", "skipped", "organ masses or metabolic rate missing")
    residual = max(mass - sum(state.get(i) for i in RMR_RATES), 0.0)
    total = sum(RMR_RATES[i] * state.get(i) for i in RMR_RATES) + RMR_RESIDUAL_RATE * residual
    stored = state.get("resting_metabolic_rate")
    diff = total / stored - 1.0
    details = {"organ_sum_kcal_day": round(total), "stored_kcal_day": round(stored), "difference": round(diff, 3),
               "basis": basis}
    msg = (f"Organ-sum metabolic rate {total:.0f} kcal/day vs ICRP basal rate {stored:.0f} ({diff:+.1%}), "
           f"residual mass from {basis}.")
    return CheckResult("rmr_organ_crosscheck", "warn" if abs(diff) > warn else "pass", msg, details)


def mifflin_crosscheck(state: State, warn: float = 0.10) -> CheckResult:
    """ICRP basal metabolic rate against the Mifflin-St Jeor equation at the ICRP reference height and mass."""
    ctx = state.context
    if ctx.age != "adult" and not (isinstance(ctx.age, (int, float)) and ctx.age >= 20):
        return CheckResult("mifflin_crosscheck", "skipped", "equation is for adults")
    h, m, bmr = state.get("reference_height"), state.get("reference_body_mass"), state.get("resting_metabolic_rate")
    if None in (h, m, bmr) or state.context.edits:
        return CheckResult("mifflin_crosscheck", "skipped", "reference height, mass or metabolic rate missing, or edited")
    pred = 10 * m + 6.25 * h - 5 * 30 + (5 if ctx.sex == "male" else -161)
    diff = bmr / pred - 1.0
    return CheckResult("mifflin_crosscheck", "warn" if abs(diff) > warn else "pass",
                       f"ICRP basal rate {bmr:.0f} kcal/day vs Mifflin-St Jeor {pred:.0f} at the reference person, "
                       f"age 30 ({diff:+.1%}).", {"icrp": round(bmr), "mifflin": round(pred)})


def icrp_red_cell_crosscheck(state: State, warn: float = 0.10) -> CheckResult:
    """Red cell volume from blood volume x haematocrit against ICRP 89 Table 2.12 (2.3 L male, 1.5 L female).
    The haematocrit comes from a different source than the blood volume, so this is a real check."""
    ref = {"male": 2.3, "female": 1.5}[state.context.sex]
    adult = state.ranges is None
    v = state.get("red_cell_volume")
    if not adult or v is None or state.context.edits:
        return CheckResult("icrp_red_cell_crosscheck", "skipped", "adult, unedited context with red cell volume needed")
    diff = v / ref - 1.0
    return CheckResult("icrp_red_cell_crosscheck", "warn" if abs(diff) > warn else "pass",
                       f"Red cell volume {v:.2f} L (blood volume x haematocrit) vs ICRP 89 Table 2.12 {ref} L ({diff:+.1%}).",
                       {"model": round(v, 3), "icrp": ref})


def solver_converged(state: State) -> CheckResult:
    if not state.solver:
        return CheckResult("solver_converged", "skipped", "solver not run (use --solve)")
    if not state.solver["converged"]:
        return CheckResult("solver_converged", "fail",
                           f"Solver did not converge in {state.solver['iterations']} iterations; "
                           f"a feedback loop has gain of 1 or more. Do not trust these values.")
    return CheckResult("solver_converged", "pass", f"Solver converged in {state.solver['iterations']} iteration(s).")


def coupling_checks(registry: Registry, couplings, sources: dict[str, dict] | None = None) -> list[CheckResult]:
    problems: list[str] = []
    for c in couplings:
        if c.target not in registry:
            problems.append(f"{c.id}: unknown target {c.target!r}")
        elif registry[c.target].kind == "derived":
            problems.append(f"{c.id}: targets derived parameter {c.target!r}")
        for pid in list(c.args.get("inputs", {})) + list(c.args.get("rates", {})) + list(c.args.get("tissues", [])) \
                + [c.args[k] for k in ("demand", "organ", "total") if k in c.args]:
            if pid not in registry:
                problems.append(f"{c.id}: unknown input {pid!r}")
        if c.source_ref and sources is not None and c.source_ref not in sources:
            problems.append(f"{c.id}: unknown source {c.source_ref!r}")
        if c.status == "referenced" and not c.source_ref:
            problems.append(f"{c.id}: status referenced but no source_ref")
    targets = [c.target for c in couplings]
    for t in sorted({t for t in targets if targets.count(t) > 1}):
        problems.append(f"{t}: more than one coupling targets it")
    return [CheckResult("coupling_data", "fail" if problems else "pass",
                        "; ".join(problems) if problems else f"{len(list(couplings))} couplings are consistent.",
                        {"problems": problems})]


def growth_continuity(registry: Registry, sex: str = "male") -> list[CheckResult]:
    """Compare the last CDC growth value (just under 20 years) with the adult baseline. A large gap is
    expected where the adult data are a different population (ANSUR II is a military sample)."""
    from .growth import load_growth
    growth = load_growth()
    out = []
    for pid in growth.parameters:
        pt = growth.at(pid, sex, 239.9)
        adult = registry[pid].baseline(sex, "adult") if pid in registry else None
        if pt is None or adult is None:
            continue
        diff = pt.value / adult - 1.0
        out.append(CheckResult(
            f"growth_continuity:{pid}:{sex}", "warn" if abs(diff) > 0.05 else "pass",
            f"CDC median at 20 y {pt.value:.1f} vs adult baseline {adult:.1f} ({diff:+.1%}). "
            f"Different populations: CDC is the general US public, ANSUR II is US Army personnel.",
            {"cdc": round(pt.value, 2), "adult": adult, "difference": round(diff, 3)}))
    return out


def bmi_consistency(state: State, tol: float = 0.08) -> CheckResult:
    """At a pediatric age the BMI computed from median height and weight should be close to the CDC
    median BMI (not identical: the median of a ratio is not the ratio of medians)."""
    from .growth import load_growth
    ctx = state.context
    if state.ranges is None or ctx.edits:
        return CheckResult("bmi_consistency", "skipped", "not an unedited pediatric context")
    ref = load_growth().check_point("body_mass_index", ctx.sex, float(ctx.age) * 12)
    mine = state.get("body_mass_index")
    if ref is None or mine is None:
        return CheckResult("bmi_consistency", "skipped", "no CDC BMI median at this age")
    diff = mine / ref - 1.0
    return CheckResult("bmi_consistency", "warn" if abs(diff) > tol else "pass",
                       f"BMI from median height and weight {mine:.2f} vs CDC median BMI {ref:.2f} ({diff:+.1%}).",
                       {"derived": round(mine, 2), "cdc": round(ref, 2)})


def state_checks(state: State) -> list[CheckResult]:
    return [tissue_mass_budget(state), segment_stature_sum(state), range_flags(state),
            sex_rule_coverage(state), rmr_organ_crosscheck(state), mifflin_crosscheck(state), icrp_red_cell_crosscheck(state), solver_converged(state),
            bmi_consistency(state)]


def modifier_checks(registry: Registry, modifiers, sources: dict[str, dict] | None = None) -> list[CheckResult]:
    """Data checks for modifiers: valid targets, known sources, and no duplicate sex rules."""
    problems: list[str] = []
    seen: dict[tuple, list[str]] = {}
    for m in modifiers:
        for pid in m.effects:
            if pid not in registry:
                problems.append(f"{m.id}: unknown parameter {pid!r}")
            elif registry[pid].kind == "derived":
                problems.append(f"{m.id}: targets derived parameter {pid!r}")
            if m.stage == "sex":
                seen.setdefault((pid, tuple(sorted(m.when.items()))), []).append(m.id)
        if m.source_ref and sources is not None and m.source_ref not in sources:
            problems.append(f"{m.id}: unknown source {m.source_ref!r}")
        if m.status == "referenced" and not m.source_ref:
            problems.append(f"{m.id}: status referenced but no source_ref")
    for m in modifiers:
        if m.stage != "sex" or m.op != "ratio" or dict(m.when) != {"sex": "female"}:
            continue
        for pid, x in m.effects.items():
            if pid in registry:
                a, b = registry[pid].baseline("male", "adult"), registry[pid].baseline("female", "adult")
                if a and b and abs(b / a - x) > 1e-3:
                    problems.append(f"{m.id}: rule ratio {x} disagrees with stored pair ratio {b / a:.5f}; "
                                    f"regenerate with tools/build_sex_rules.py")
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    for (pid, _), ids in dups.items():
        problems.append(f"duplicate sex rules for {pid!r}: {', '.join(ids)} (would apply twice)")
    return [CheckResult("modifier_data", "fail" if problems else "pass",
                        "; ".join(problems) if problems else f"{len(list(modifiers))} modifiers are consistent.",
                        {"problems": problems})]


def dataset_checks(registry: Registry, sources: dict[str, dict], extra_used=()) -> list[CheckResult]:
    """extra_used: source ids cited outside the registry (modifiers, couplings, validation files)."""
    results: list[CheckResult] = []

    bad_refs = []
    for p in registry:
        for ref in (p.source_ref, p.range_source_ref):
            if ref and ref not in sources:
                bad_refs.append(f"{p.id}->{ref}")
    for sid, rec in sources.items():
        primary = rec.get("primary")
        if primary and primary not in sources:
            bad_refs.append(f"{sid}->primary:{primary}")
    results.append(CheckResult(
        "source_refs_exist", "fail" if bad_refs else "pass",
        f"Unknown source ids: {', '.join(bad_refs)}" if bad_refs else "Every source_ref exists.",
        {"bad": bad_refs}))

    no_ref = [p.id for p in registry if p.status == "referenced" and not p.source_ref]
    results.append(CheckResult(
        "referenced_have_source", "fail" if no_ref else "pass",
        f"Status 'referenced' without source_ref: {', '.join(no_ref)}" if no_ref
        else "Every referenced value names a source.", {"bad": no_ref}))

    used = {r for p in registry for r in (p.source_ref, p.range_source_ref) if r} | set(extra_used)
    used |= {sources[r]["primary"] for r in used if r in sources and sources[r].get("primary")}
    unused = sorted(set(sources) - used)
    results.append(CheckResult(
        "unused_sources", "warn" if unused else "pass",
        f"Sources not used by any parameter: {', '.join(unused)}" if unused
        else "Every source is used.", {"unused": unused}))
    return results
