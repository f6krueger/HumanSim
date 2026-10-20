"""Couplings: the edges of the parameter network.

A coupling says how one stored parameter (the target) responds when others change. Every coupling
is **reference-calibrated**: it returns a ratio g, and the solved value is

    target = reference_target * g

where the reference state is the same context with no user edits. So the unedited model always
reproduces its stored values exactly, and a coupling only describes how a value *moves*.

Each coupling is labelled:

  kind   mechanistic  a physical or physiological law (Fick principle, metabolic sum)
         empirical    a fitted relationship
         assumed      a modelling assumption with no source; flagged on every value it moves
  role   accounting   bookkeeping that must hold (masses add up, metabolic rate is a sum);
                      applies at once, with no deadband
         compensation the body adjusting to a change; held back by a deadband so that changes
                      within normal variation leave the target alone
  deadband   a float (fraction either side of the reference) or "ref_range" (each driving input's
             own 2.5-97.5 percentile range for the context sex). Inside the band the driving ratio
             counts as 1; beyond it the response starts from the band edge, so it is continuous.
  validity   [low, high] range of the (undamped) driving ratio the relationship is believed to
             hold for. Outside it the solved value is flagged outside_coupling_validity.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping

import yaml

from . import derived as _derived
from .registry import DATA_DIR, Registry, RegistryError

KINDS = ("mechanistic", "empirical", "assumed")
ROLES = ("accounting", "compensation")


@dataclass(frozen=True)
class Coupling:
    id: str
    target: str
    fn: str
    args: Mapping = field(default_factory=dict)
    kind: str = "assumed"
    role: str = "compensation"
    deadband: object = 0.0  # float | "ref_range"
    validity: tuple[float, float] | None = None
    status: str = "unsourced"
    source_ref: str | None = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise RegistryError(f"Coupling {self.id!r}: kind must be one of {KINDS}")
        if self.role not in ROLES:
            raise RegistryError(f"Coupling {self.id!r}: role must be one of {ROLES}")
        if self.fn not in FUNCTIONS:
            raise RegistryError(f"Coupling {self.id!r}: unknown function {self.fn!r}")
        if self.role == "accounting" and self.deadband:
            raise RegistryError(f"Coupling {self.id!r}: accounting couplings have no deadband")
        if self.deadband != "ref_range" and not isinstance(self.deadband, (int, float)):
            raise RegistryError(f"Coupling {self.id!r}: deadband must be a number or 'ref_range'")


class SolveCtx:
    """What a coupling function sees: current and reference values, and the deadband rule."""

    def __init__(self, registry: Registry, sex: str, cur: dict, ref: dict, pinned: set[str],
                 coupling: Coupling):
        self.registry, self.sex, self._cur, self._ref = registry, sex, cur, ref
        self.pinned, self.coupling = pinned, coupling

    def cur(self, pid: str):
        return self._cur.get(pid)

    def ref(self, pid: str):
        return self._ref.get(pid)

    def ratio(self, pid: str) -> float | None:
        c, r = self.cur(pid), self.ref(pid)
        if c is None or r in (None, 0):
            return None
        return c / r

    def dead(self, ratio: float, pid: str | None = None) -> float:
        """Apply this coupling's deadband to a driving ratio."""
        lo, hi = self.band(pid)
        if lo <= ratio <= hi:
            return 1.0
        return ratio / hi if ratio > hi else ratio / lo

    def band(self, pid: str | None) -> tuple[float, float]:
        d = self.coupling.deadband
        if d == "ref_range":
            rng = self.registry[pid].ref_range.get(self.sex) if pid else None
            base = self.ref(pid) if pid else None
            if not rng or not base:
                return (1.0, 1.0)  # range undefined: no deadband, never a guess at one
            return (rng[0] / base, rng[1] / base)
        return (1.0 - float(d), 1.0 + float(d))


# Each function returns (g, raw) or None when an input is missing:
#   g    the ratio to multiply the reference target by (after deadband where the role allows it)
#   raw  the undamped driving ratio, used for the validity check
def power_ratio(ctx: SolveCtx, args: Mapping):
    g = raw = 1.0
    for pid, exp in args["inputs"].items():
        r = ctx.ratio(pid)
        if r is None:
            return None
        raw *= r ** exp
        g *= ctx.dead(r, pid) ** exp
    return g, raw


def _cao2(ctx: SolveCtx, which: str) -> float | None:
    src = ctx._cur if which == "cur" else ctx._ref
    need = {k: src.get(k) for k in ("hemoglobin", "sao2", "pao2")}
    if any(v is None for v in need.values()):
        return None
    return _derived.compute(need)["arterial_o2_content"]


def fick_cardiac_output(ctx: SolveCtx, args: Mapping):
    """Fick principle: CO = O2 consumption / arterial-venous O2 difference. With the extraction held
    constant, CO scales with metabolic rate and inversely with arterial O2 content."""
    r = ctx.ratio(args["demand"])
    cur, ref = _cao2(ctx, "cur"), _cao2(ctx, "ref")
    if r is None or not cur or not ref:
        return None
    raw = r * ref / cur
    return ctx.dead(raw), raw


def _rmr_parts(values: Mapping, args: Mapping):
    total = values.get(args["total"])
    rates = args["rates"]
    masses = {pid: values.get(pid) for pid in rates}
    if total is None or any(m is None for m in masses.values()):
        return None
    residual = total - sum(masses.values())
    return {pid: rates[pid] * m for pid, m in masses.items()}, args["residual_rate"] * max(residual, 0.0)


def rmr_organ_sum(ctx: SolveCtx, args: Mapping):
    """Resting metabolic rate as the sum of tissue mass times its specific metabolic rate (Wang 2010)."""
    cur, ref = _rmr_parts(ctx._cur, args), _rmr_parts(ctx._ref, args)
    if cur is None or ref is None:
        return None
    s_cur, s_ref = sum(cur[0].values()) + cur[1], sum(ref[0].values()) + ref[1]
    return s_cur / s_ref, s_cur / s_ref


def metabolic_share(ctx: SolveCtx, args: Mapping):
    """An organ's blood flow follows its share of oxygen use (extraction assumed constant)."""
    cur, ref = _rmr_parts(ctx._cur, args), _rmr_parts(ctx._ref, args)
    if cur is None or ref is None:
        return None
    organ = args["organ"]
    share_cur = cur[0][organ] / (sum(cur[0].values()) + cur[1])
    share_ref = ref[0][organ] / (sum(ref[0].values()) + ref[1])
    raw = share_cur / share_ref
    return ctx.dead(raw), raw


def body_mass_balance(ctx: SolveCtx, args: Mapping):
    """Body mass = reference mass scaled with height (BMI-preserving exponent) + the change in any
    tissue the user edited directly. Only edited tissues are added: tissues that move because of
    other couplings are already part of the mass those couplings respond to."""
    h = ctx.ratio("height")
    ref_mass = ctx.ref(ctx.coupling.target)
    if h is None or not ref_mass:
        return None
    mass = ref_mass * ctx.dead(h, "height") ** args["height_exponent"]
    for pid in args.get("tissues", ()):
        if pid in ctx.pinned and ctx.cur(pid) is not None and ctx.ref(pid) is not None:
            mass += ctx.cur(pid) - ctx.ref(pid)
    return mass / ref_mass, mass / ref_mass


FUNCTIONS: dict[str, Callable] = {
    "power_ratio": power_ratio,
    "fick_cardiac_output": fick_cardiac_output,
    "rmr_organ_sum": rmr_organ_sum,
    "metabolic_share": metabolic_share,
    "body_mass_balance": body_mass_balance,
}


def load_couplings(registry: Registry | None = None, data_dir: str | Path | None = None) -> list[Coupling]:
    """Load data/couplings/*.yaml. A record with `targets` expands to one coupling per target, and a
    `select` filter ({group, unit, exclude}) picks targets from the registry."""
    root = Path(data_dir) if data_dir else DATA_DIR
    out: list[Coupling] = []
    for path in sorted((root / "couplings").glob("*.yaml")):
        doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for raw in doc.get("couplings") or []:
            raw = dict(raw)
            targets = raw.pop("targets", None)
            select = raw.pop("select", None)
            if raw.get("validity") is not None:
                raw["validity"] = tuple(raw["validity"])
            if select:
                if registry is None:
                    raise RegistryError(f"{path.name}: a 'select' needs the registry")
                targets = [p.id for p in registry
                           if p.kind != "derived" and p.group == select.get("group", p.group)
                           and p.unit == select.get("unit", p.unit)
                           and p.id not in select.get("exclude", ())]
            if targets is None:
                out.append(Coupling(**raw))
                continue
            base = raw.pop("id")
            raw.pop("target", None)
            for t in targets:
                out.append(Coupling(id=f"{base}:{t}", target=t, **raw))
    ids = [c.id for c in out]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise RegistryError(f"Duplicate coupling ids: {', '.join(dup)}")
    return out
