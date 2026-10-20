"""A small generic ODE engine (classical fourth-order Runge-Kutta, fixed step).

Models are plain functions f(t, y, u) -> dy/dt over dicts of named variables, so a model is a few
lines and the engine knows nothing about physiology. Time is in days throughout HumanSim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Mapping

Vector = dict[str, float]
Rhs = Callable[[float, Vector, float], Vector]


@dataclass
class Trajectory:
    t: list[float]
    y: list[Vector]
    converged: bool = True  # False if a value became non-finite

    def series(self, name: str) -> list[float]:
        return [s[name] for s in self.y]

    def final(self) -> Vector:
        return self.y[-1]


def _add(y: Mapping[str, float], k: Mapping[str, float], h: float) -> Vector:
    return {n: y[n] + h * k[n] for n in y}


def rk4(f: Rhs, y0: Mapping[str, float], t_end: float, dt: float, t0: float = 0.0,
        u: Callable[[float], float] | float = 0.0) -> Trajectory:
    """Integrate dy/dt = f(t, y, u(t)) from t0 to t_end with step dt. `u` is an input signal
    (a function of time or a constant)."""
    if dt <= 0 or t_end < t0:
        raise ValueError("need dt > 0 and t_end >= t0")
    uf = u if callable(u) else (lambda _t, _v=float(u): _v)
    n = max(1, round((t_end - t0) / dt))
    h = (t_end - t0) / n
    t, y = t0, {k: float(v) for k, v in y0.items()}
    ts, ys = [t], [dict(y)]
    for _ in range(n):
        try:
            k1 = f(t, y, uf(t))
            k2 = f(t + h / 2, _add(y, k1, h / 2), uf(t + h / 2))
            k3 = f(t + h / 2, _add(y, k2, h / 2), uf(t + h / 2))
            k4 = f(t + h, _add(y, k3, h), uf(t + h))
        except (OverflowError, ZeroDivisionError):
            return Trajectory(ts, ys, converged=False)  # the model ran away
        y = {name: y[name] + h / 6 * (k1[name] + 2 * k2[name] + 2 * k3[name] + k4[name]) for name in y}
        t += h
        ts.append(t)
        ys.append(dict(y))
        if not all(math.isfinite(v) for v in y.values()):
            return Trajectory(ts, ys, converged=False)
    return Trajectory(ts, ys)
