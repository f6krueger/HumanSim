"""Validation targets: published statistics that derived values must reproduce.

These files are not baseline data. The Phase 2 sex rules derive the female body from the male
baseline, and the female statistics here are what that derivation is checked against.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from .registry import DATA_DIR, Registry


def load_validation(name: str, data_dir: str | Path | None = None) -> dict:
    root = Path(data_dir) if data_dir else DATA_DIR
    path = root / "validation" / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"No validation file {path}")
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


# Convention: a measurement is "scaling-only" if the female/male ratio is within 3% of what
# geometric similarity predicts from the stature ratio. This threshold is a modelling choice.
SCALING_TOLERANCE = 0.03


def scaling_expectation(unit: str, height_ratio: float) -> float | None:
    """Female/male ratio predicted by size alone (geometric similarity), or None if not defined."""
    if unit in ("cm", "mm"):
        return height_ratio  # linear dimensions scale with stature
    if unit == "kg":
        return height_ratio ** 3  # mass scales with volume
    return None


def classify_ratio(ratio: float, expected: float | None, tol: float = SCALING_TOLERANCE) -> str:
    """'scaling' if size alone explains the ratio, otherwise 'anatomical'."""
    if expected is None:
        return "anatomical"
    return "scaling" if abs(ratio / expected - 1.0) <= tol else "anatomical"


def scaling_analysis(doc: dict, registry: Registry) -> list[dict]:
    """For each validation parameter: ratio, size-only prediction, deviation, classification."""
    params = doc["parameters"]
    height_ratio = params["height"]["mean"] / registry["height"].baseline("male", "adult")
    out = []
    for pid, rec in params.items():
        male = registry[pid].baseline("male", "adult")
        ratio = rec["mean"] / male
        expected = scaling_expectation(rec["unit"], height_ratio)
        out.append({
            "id": pid, "unit": rec["unit"], "ratio": ratio, "expected": expected,
            "deviation": None if expected is None else ratio / expected - 1.0,
            "type": classify_ratio(ratio, expected),
        })
    return out


def validation_problems(doc: dict, registry: Registry) -> list[str]:
    """Return a list of problems (empty when the file is consistent with the registry)."""
    problems = []
    for pid, rec in (doc.get("parameters") or {}).items():
        if pid not in registry:
            problems.append(f"{pid}: not a registry parameter")
            continue
        if rec.get("unit") != registry[pid].unit:
            problems.append(f"{pid}: unit {rec.get('unit')!r} != registry {registry[pid].unit!r}")
        if not (rec["p2_5"] < rec["mean"] < rec["p97_5"]):
            problems.append(f"{pid}: mean outside its own 2.5-97.5 range")
        if rec["sd"] <= 0 or rec["n"] <= 0:
            problems.append(f"{pid}: non-positive sd or n")
    return problems
