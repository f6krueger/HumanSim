"""Consistency checks.

state_checks() inspects one resolved State. dataset_checks() inspects the registry
against the sources sidecar. A check with too little data reports "skipped", never "pass".
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .context import State
from .registry import Registry

BLOOD_DENSITY_KG_PER_L = 1.06  # approximate
MAX_UNASSIGNED_FRACTION = 0.10  # above this, the tissue budget warns

# Tissues whose masses are expected to add up to (most of) body mass.
TISSUE_MASSES = (
    "skeletal_muscle_mass", "skeleton_mass", "adipose_tissue_mass", "skin_mass",
    "liver_mass", "brain_mass", "heart_mass", "kidneys_mass", "lungs_mass",
    "spleen_mass", "pancreas_mass", "thyroid_mass",
)


@dataclass
class CheckResult:
    id: str
    status: str  # pass | warn | fail | skipped
    message: str
    details: dict = field(default_factory=dict)


def tissue_mass_budget(state: State) -> CheckResult:
    mass = state.get("body_mass")
    if mass is None:
        return CheckResult("tissue_mass_budget", "skipped", "body_mass has no value")
    missing = [i for i in TISSUE_MASSES if state.get(i) is None]
    total = sum(state.get(i) or 0.0 for i in TISSUE_MASSES)
    blood = state.get("blood_volume")
    if blood is not None:
        total += blood * BLOOD_DENSITY_KG_PER_L
    else:
        missing.append("blood_volume")
    remainder = mass - total
    details = {"listed_total_kg": round(total, 3), "unassigned_kg": round(remainder, 3),
               "missing": missing}
    if total > mass:
        return CheckResult("tissue_mass_budget", "fail",
                           f"Listed tissues sum to {total:.2f} kg, more than body mass {mass:.2f} kg.",
                           details)
    if missing:
        return CheckResult("tissue_mass_budget", "warn",
                           f"Partial budget: {len(missing)} tissue(s) have no value; "
                           f"{remainder:.2f} kg unassigned.", details)
    if remainder / mass > MAX_UNASSIGNED_FRACTION:
        return CheckResult("tissue_mass_budget", "warn",
                           f"{remainder:.1f} kg ({remainder / mass:.0%}) of body mass is not covered by "
                           f"the listed tissues. Most likely the organ and tissue masses (ICRP reference "
                           f"man, 73 kg) and body mass ({mass:.1f} kg) come from different populations; "
                           f"it does not mean anatomy is missing.", details)
    return CheckResult("tissue_mass_budget", "pass",
                       f"Listed tissues sum to {total:.2f} kg; {remainder:.2f} kg unassigned "
                       f"(gut contents, bone marrow detail, double-counted blood in organs, etc.).",
                       details)


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


def female_rule_coverage(state: State) -> CheckResult:
    if state.context.sex != "female":
        return CheckResult("female_rule_coverage", "skipped", "not a female context")
    no_rule = sorted(i for i, f in state.flags.items() if "no_sex_rule" in f)
    if no_rule:
        shown = ", ".join(no_rule[:8]) + (f", and {len(no_rule) - 8} more" if len(no_rule) > 8 else "")
        return CheckResult(
            "female_rule_coverage", "warn",
            f"{len(no_rule)} parameter(s) have a male value but no sex rule, so they are "
            f"empty (or flagged, with keep_unadjusted): {shown}", {"no_sex_rule": no_rule})
    return CheckResult("female_rule_coverage", "pass", "Every parameter with a male value has a sex rule.")


def state_checks(state: State) -> list[CheckResult]:
    return [tissue_mass_budget(state), segment_stature_sum(state), range_flags(state),
            female_rule_coverage(state)]


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
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    for (pid, _), ids in dups.items():
        problems.append(f"duplicate sex rules for {pid!r}: {', '.join(ids)} (would apply twice)")
    return [CheckResult("modifier_data", "fail" if problems else "pass",
                        "; ".join(problems) if problems else f"{len(list(modifiers))} modifiers are consistent.",
                        {"problems": problems})]


def dataset_checks(registry: Registry, sources: dict[str, dict]) -> list[CheckResult]:
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

    used = {r for p in registry for r in (p.source_ref, p.range_source_ref) if r}
    used |= {sources[r]["primary"] for r in used if r in sources and sources[r].get("primary")}
    unused = sorted(set(sources) - used)
    results.append(CheckResult(
        "unused_sources", "warn" if unused else "pass",
        f"Sources not used by any parameter: {', '.join(unused)}" if unused
        else "Every source is used.", {"unused": unused}))
    return results
