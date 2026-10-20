"""Pharmacokinetics (Phase 6): what a dose does to the concentration of a drug in this body.

Scope, plainly: a **one-compartment model with intravenous dosing**, not a full PBPK model. Each
drug has a reference clearance and volume of distribution per kg from a cited label. For a given
person they are scaled by the model's own values, the same reference-calibrated way as the solver:

    CL = CL_ref_kg * reference_mass * prod((driver / driver_ref) ** exponent)
    V  = V_ref_kg  * reference_mass * prod((driver / driver_ref) ** exponent)

where the reference is the same sex and age with no edits, hormone or other imposed change, and the
drivers (liver mass, hepatic blood flow, body mass) are read from the resolved or solved State. So
a larger liver, a lower cardiac output or a heavier body changes the curve. The exponents are
assumptions and are reported. Time is in HOURS here (the hormone models use days).

Not modelled: oral absorption (no sourced bioavailability or absorption rate), renal drugs (no
sourced kidney clearance values yet), metabolites, protein binding, drug effects on the body, and
variability between people.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import yaml

from .context import Context, State, resolve
from .registry import DATA_DIR, Registry, RegistryError

DRIVERS = ("liver_mass", "hepatic_blood_flow", "body_mass")


@dataclass(frozen=True)
class Drug:
    id: str
    name: str
    routes: tuple[str, ...]
    clearance_l_h_kg: float
    volume_l_kg: float
    label_half_life_h: float | None
    status: str
    source_ref: str | None
    clearance_scaling: Mapping[str, float]
    volume_scaling: Mapping[str, float]
    scaling_note: str = ""
    genotype_effects: tuple = ()  # ({gene, auc_ratio_by_activity_score, source_ref, note}, ...)

    def __post_init__(self) -> None:
        if self.clearance_l_h_kg <= 0 or self.volume_l_kg <= 0:
            raise RegistryError(f"Drug {self.id!r}: clearance and volume must be positive")
        for scaling in (self.clearance_scaling, self.volume_scaling):
            bad = set(scaling) - set(DRIVERS)
            if bad:
                raise RegistryError(f"Drug {self.id!r}: unknown scaling driver(s) {sorted(bad)}")
        if self.status == "referenced" and not self.source_ref:
            raise RegistryError(f"Drug {self.id!r}: referenced values need a source")


def load_drugs(data_dir: str | Path | None = None) -> dict[str, Drug]:
    root = Path(data_dir) if data_dir else DATA_DIR
    out: dict[str, Drug] = {}
    for path in sorted((root / "drugs").glob("*.yaml")):
        for raw in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("drugs") or []:
            raw = dict(raw)
            d = Drug(id=raw["id"], name=raw["name"], routes=tuple(raw["routes"]),
                     clearance_l_h_kg=float(raw["clearance_l_h_kg"]), volume_l_kg=float(raw["volume_l_kg"]),
                     label_half_life_h=raw.get("label_half_life_h"), status=raw.get("status", "unsourced"),
                     source_ref=raw.get("source_ref"), clearance_scaling=dict(raw.get("clearance_scaling") or {}),
                     volume_scaling=dict(raw.get("volume_scaling") or {}), scaling_note=raw.get("scaling_note", ""),
                     genotype_effects=tuple(raw.get("genotype_effects") or ()))
            if d.id in out:
                raise RegistryError(f"Duplicate drug id {d.id!r}")
            out[d.id] = d
    return out


@dataclass(frozen=True)
class Dose:
    amount_mg: float
    time_h: float = 0.0
    duration_h: float = 0.0  # 0 = bolus (given over a minute, see below); >0 = constant-rate infusion


@dataclass
class PkResult:
    drug: Drug
    clearance_l_h: float
    volume_l: float
    times_h: list[float]
    conc_mg_l: list[float]
    flags: set[str] = field(default_factory=set)
    notes: list[str] = field(default_factory=list)

    @property
    def half_life_h(self) -> float:
        return math.log(2) * self.volume_l / self.clearance_l_h

    @property
    def cmax_mg_l(self) -> float:
        return max(self.conc_mg_l)

    def auc_mg_h_l(self) -> float:
        """Area under the curve over the simulated span (trapezoid)."""
        t, c = self.times_h, self.conc_mg_l
        return sum((t[i + 1] - t[i]) * (c[i] + c[i + 1]) / 2 for i in range(len(t) - 1))


def _driver(state: State, name: str) -> float:
    v = state.values.get(name)
    if v is None:
        raise ValueError(f"{name} has no value in this state, so the drug cannot be scaled")
    return float(v)


def _scale(state: State, ref: State, scaling: Mapping[str, float]) -> float:
    f = 1.0
    for name, exp in scaling.items():
        f *= (_driver(state, name) / _driver(ref, name)) ** exp
    return f


def genotype_clearance_factor(drug: Drug, genotype: Mapping[str, str], genes=None) -> tuple[float, list[str], set[str]]:
    """(factor on clearance, notes, flags) from the person's genotype. 1.0 for the reference (*1/*1).
    A gene the person has a genotype for but the drug has no data on changes nothing and is noted."""
    from .genetics import exposure_ratio, load_genes

    genes = genes or load_genes()
    factor, notes, flags = 1.0, [], set()
    for gid, dip in genotype.items():
        if gid not in genes:
            raise KeyError(f"Unknown gene {gid!r}; known: {', '.join(genes)}")
        score = genes[gid].activity_score(dip)
        effects = [e for e in drug.genotype_effects if e["gene"] == gid]
        if not effects:
            notes.append(f"{gid} {dip}: no data on {drug.id}, so no effect applied")
            continue
        for e in effects:
            ratio, interpolated = exposure_ratio(e["auc_ratio_by_activity_score"], score)
            factor /= ratio
            flags |= {"genotype_effect"} | ({"interpolated"} if interpolated else set())
            notes.append(f"{gid} {dip} (activity score {score:g}, {genes[gid].phenotype(dip)}): "
                         f"AUC x{ratio:.3g}, clearance /{ratio:.3g}" + (" (interpolated)" if interpolated else ""))
    return factor, notes, flags


def pk_parameters(registry: Registry, state: State, drug: Drug, modifiers=()) -> tuple[float, float, float]:
    """(clearance L/h, volume L, reference mass kg) for this state."""
    ctx = state.context
    ref = resolve(registry, Context(sex=ctx.sex, age=ctx.age, cycle_phase=ctx.cycle_phase,
                                    growth_z=ctx.growth_z, keep_unadjusted=ctx.keep_unadjusted), modifiers)
    ref_mass = ref.values.get("reference_body_mass") or ref.values.get("body_mass")
    if ref_mass is None:
        raise ValueError("no body mass for this context")
    cl = drug.clearance_l_h_kg * ref_mass * _scale(state, ref, drug.clearance_scaling)
    cl *= genotype_clearance_factor(drug, ctx.genotype)[0]
    v = drug.volume_l_kg * ref_mass * _scale(state, ref, drug.volume_scaling)
    return cl, v, ref_mass


def simulate(registry: Registry, state: State, drug: Drug, doses: list[Dose], hours: float,
             modifiers=(), dt: float = 0.01) -> PkResult:
    """Concentration of the drug in plasma (mg/L) over `hours`, for intravenous doses."""
    if "iv" not in drug.routes:
        raise ValueError(f"{drug.id}: no intravenous route in the data")
    if not doses or hours <= 0:
        raise ValueError("need at least one dose and a positive duration")
    cl, v, _ = pk_parameters(registry, state, drug, modifiers)
    k = cl / v
    bolus_window = 1.0 / 60.0  # a bolus is delivered over one minute
    spans = [(d.time_h, d.duration_h if d.duration_h > 0 else bolus_window, d.amount_mg) for d in doses]

    def conc(t: float) -> float:
        """Exact superposition of constant-rate infusions into one compartment (no numerical error)."""
        total = 0.0
        for t0, dur, amount in spans:
            if t <= t0:
                continue
            r = amount / dur
            if t <= t0 + dur:
                total += r / cl * (1 - math.exp(-k * (t - t0)))
            else:
                total += r / cl * (1 - math.exp(-k * dur)) * math.exp(-k * (t - t0 - dur))
        return total

    grid = {round(i * dt, 9) for i in range(int(hours / dt) + 1)} | {hours}
    for t0, dur, _ in spans:
        grid |= {t for t in (t0, t0 + dur) if 0 <= t <= hours}
    times = sorted(grid)
    res = PkResult(drug, cl, v, times, [conc(t) for t in times])
    res.flags = {"assumed_scaling", "one_compartment"} | ({"unsourced"} if drug.status != "referenced" else set())
    if any(f in state.flags.get(p, set()) for p in drug.clearance_scaling for f in ("edited", "modified", "solved")):
        res.flags.add("scaled_from_changed_body")
    res.notes.append(drug.scaling_note)
    _, gnotes, gflags = genotype_clearance_factor(drug, state.context.genotype)
    res.notes.extend(gnotes)
    res.flags |= gflags
    return res


def half_life_check(registry: Registry, drug: Drug, state: State, modifiers=()) -> tuple[float, float | None]:
    """(model half-life, label half-life) for the unedited person, so the gap is visible."""
    cl, v, _ = pk_parameters(registry, state, drug, modifiers)
    return math.log(2) * v / cl, drug.label_half_life_h
