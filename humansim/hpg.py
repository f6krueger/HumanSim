"""Male hypothalamic-pituitary-gonadal (HPG) feedback model.

Three variables, in registry units: LH (mIU/mL), testicular responsiveness C (dimensionless, 1 =
reference) and total testosterone T (ng/dL).

    dLH/dt = (LH* / tau_LH) * ((T / T*) ** -n  -  LH / LH*)     LH is driven up when T falls
    dC/dt  = (LH / LH*  -  C) / tau_C                           the testes follow LH over weeks
    dT/dt  = (T* * C  -  T) / tau_T  +  u(t)                    production minus clearance, plus dosing

LH* and T* are the stored male adult reference values, so with no input the model stays exactly at
the reference. Everything else here is an ASSUMPTION, not data: the structure (a simple negative
feedback loop with a slow gonadal term), the feedback strength n and the three time constants.
They are set so the model behaves plausibly (hours for hormone clearance, weeks for the testes), and
every result carries `assumed: True`. Male only: the female model is the cycle profile (cycle.py).
Hormone-effect couplings come from sourced data elsewhere; this module does not change them.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dynamics import Trajectory, rk4
from .registry import Registry

ASSUMED = {
    "feedback_strength_n": 2.0,   # how hard T suppresses LH (assumed)
    "tau_lh_days": 0.05,          # LH clearance, about an hour (assumed)
    "tau_t_days": 0.05,           # T clearance, about an hour (assumed)
    "tau_c_days": 30.0,           # testicular responsiveness follows LH over weeks (assumed)
}


@dataclass
class HpgResult:
    trajectory: Trajectory
    reference: dict[str, float]
    assumed: bool = True

    def final(self) -> dict[str, float]:
        return self.trajectory.final()


def reference_levels(registry: Registry) -> tuple[float, float]:
    lh = registry["lh"].baseline("male", "adult")
    t = registry["testosterone_total"].baseline("male", "adult")
    if lh is None or t is None:
        raise ValueError("male adult LH and testosterone are needed")
    return lh, t


def simulate(registry: Registry, days: float, dose_rate=0.0, start: dict[str, float] | None = None,
             dt: float = 0.01, params: dict[str, float] | None = None) -> HpgResult:
    """Run the model. `dose_rate` is exogenous testosterone input in ng/dL per day (a number or a
    function of time in days). `start` overrides the starting levels (default: the reference)."""
    p = {**ASSUMED, **(params or {})}
    lh_ref, t_ref = reference_levels(registry)
    n = p["feedback_strength_n"]

    def f(_t, y, u):
        t_level = max(y["T"], 1e-9)
        return {
            "LH": lh_ref / p["tau_lh_days"] * ((t_level / t_ref) ** -n - y["LH"] / lh_ref),
            "C": (y["LH"] / lh_ref - y["C"]) / p["tau_c_days"],
            "T": (t_ref * y["C"] - y["T"]) / p["tau_t_days"] + u,
        }

    y0 = {"LH": lh_ref, "C": 1.0, "T": t_ref}
    y0.update(start or {})
    return HpgResult(rk4(f, y0, days, dt, u=dose_rate), {"LH": lh_ref, "T": t_ref})
