"""Menstrual-cycle hormone profile for an average adult woman.

The registry stores LH, FSH, oestradiol and progesterone for follicular, midcycle and luteal phases
(single-lab interval midpoints). This module turns those phase values into a profile over the cycle
by linear interpolation between phase anchors placed at the middle of each phase (days 7, 14 and 21
of a 28-day cycle, repeating). It is a smooth sketch: the real midcycle LH peak is brief and sharp,
so the interpolated profile understates peaks and day-to-day shape. Every value is flagged
`interpolated`. A phase with no stored value for a hormone (progesterone has no midcycle entry) is
simply not an anchor, so its neighbours are joined directly.
"""

from __future__ import annotations

from .context import Context
from .registry import Registry

HORMONES = ("lh", "fsh", "estradiol", "progesterone")
ANCHOR_DAYS = {"follicular": 7.0, "midcycle": 14.0, "luteal": 21.0}
CYCLE_LENGTH = 28.0


def phase_of(day: float) -> str:
    d = (day - 1.0) % CYCLE_LENGTH + 1.0
    return "follicular" if d < 13 else "midcycle" if d < 16 else "luteal"


def _anchors(registry: Registry, pid: str) -> list[tuple[float, float]]:
    p = registry[pid]
    pts = [(day, p.baseline("female", f"adult:{ph}")) for ph, day in ANCHOR_DAYS.items()]
    return [(d, v) for d, v in pts if v is not None]


def level(registry: Registry, pid: str, day: float) -> float | None:
    """Level of one hormone on a cycle day (1-28, wraps), or None if it has no phase values."""
    pts = _anchors(registry, pid)
    if not pts:
        return None
    d = (day - 1.0) % CYCLE_LENGTH + 1.0
    # periodic: repeat the anchors one cycle before and after
    ext = [(a - CYCLE_LENGTH, v) for a, v in pts] + pts + [(a + CYCLE_LENGTH, v) for a, v in pts]
    for (a0, v0), (a1, v1) in zip(ext, ext[1:]):
        if a0 <= d <= a1:
            return v0 + (v1 - v0) * (d - a0) / (a1 - a0)
    return None


def profile(registry: Registry, days: list[float] | None = None) -> dict[float, dict[str, float]]:
    days = days if days is not None else [float(d) for d in range(1, 29)]
    out = {}
    for d in days:
        out[d] = {h: v for h in HORMONES if (v := level(registry, h, d)) is not None}
    return out


def cycle_context(registry: Registry, day: float, **kwargs) -> Context:
    """A female Context with the hormones set to the interpolated profile on `day` and the cycle
    phase set to the phase that day falls in."""
    levels = profile(registry, [day])[day]
    return Context(sex="female", hormone_state=levels, hormone_source=f"cycle profile day {day:g}",
                   cycle_phase=phase_of(day), **kwargs)
