"""Growth curves: the age axis for children and adolescents (birth to 20 years).

Values come from LMS tables (CDC 2000). For a sex and age in months, L, M and S are interpolated
linearly between table rows (the CDC's recommended procedure), and a value at z-score z is

    M * (1 + L*S*z) ** (1/L)        (M * exp(S*z) when L is 0)

z = 0 is the median, which stands for "the average child". Anything computed from these tables is
flagged `growth_curve` and `interpolated`. Parameters with no table have no value at a pediatric
age: adult values are never used for a child.
"""

from __future__ import annotations

import csv
import math
from bisect import bisect_right
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

from .registry import DATA_DIR, RegistryError

ADULT_AGE_YEARS = 20.0  # the CDC tables end at 240 months
Z_P3, Z_P97 = -1.880793608, 1.880793608  # standard normal 3rd and 97th percentiles


@dataclass(frozen=True)
class GrowthPoint:
    value: float
    sd: float   # half the distance between z=-1 and z=+1
    low: float  # 3rd percentile
    high: float  # 97th percentile
    source_ref: str


def lms_value(L: float, M: float, S: float, z: float) -> float:
    if abs(L) < 1e-9:
        return M * math.exp(S * z)
    return M * (1.0 + L * S * z) ** (1.0 / L)


class Table:
    """Rows of one LMS table for one sex, sorted by age in months."""

    def __init__(self, rows: list[tuple[float, float, float, float]]):
        self.rows = sorted(rows)
        self.ages = [r[0] for r in self.rows]

    def lms(self, months: float) -> tuple[float, float, float] | None:
        if not self.ages or months < self.ages[0] or months > self.ages[-1]:
            return None
        i = bisect_right(self.ages, months) - 1
        if i >= len(self.rows) - 1:
            return self.rows[-1][1:]
        a, b = self.rows[i], self.rows[i + 1]
        t = (months - a[0]) / (b[0] - a[0])
        return tuple(a[k] + t * (b[k] - a[k]) for k in (1, 2, 3))


def read_table(path: Path) -> dict[str, Table]:
    out: dict[int, list] = {1: [], 2: []}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out[int(r["Sex"])].append((float(r["Agemos"]), float(r["L"]), float(r["M"]), float(r["S"])))
    return {"male": Table(out[1]), "female": Table(out[2])}


class GrowthData:
    def __init__(self, curves: dict[str, dict], tables: dict[str, dict[str, Table]], checks: dict[str, dict]):
        self._curves, self._tables, self.checks = curves, tables, checks

    def has(self, pid: str) -> bool:
        return pid in self._curves

    @property
    def parameters(self) -> list[str]:
        return list(self._curves)

    def at(self, pid: str, sex: str, months: float, z: float = 0.0) -> GrowthPoint | None:
        spec = self._curves.get(pid)
        if spec is None:
            return None
        for t in spec["tables"]:
            if t["from_months"] <= months <= t["to_months"]:
                lms = self._tables[t["file"]][sex].lms(months)
                if lms is None:
                    continue
                L, M, S = lms
                return GrowthPoint(lms_value(L, M, S, z), (lms_value(L, M, S, 1) - lms_value(L, M, S, -1)) / 2,
                                   lms_value(L, M, S, Z_P3), lms_value(L, M, S, Z_P97), spec["source_ref"])
        return None

    def check_point(self, pid: str, sex: str, months: float) -> float | None:
        """Median of a check table (used to validate derived values), or None."""
        spec = self.checks.get(pid)
        if not spec or not spec["from_months"] <= months <= spec["to_months"]:
            return None
        lms = self._tables[spec["file"]][sex].lms(months)
        return None if lms is None else lms[1]


@lru_cache(maxsize=4)
def _load(root: str) -> GrowthData:
    base = Path(root) / "growth"
    path = base / "curves.yaml"
    if not path.exists():
        return GrowthData({}, {}, {})
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    curves, checks, tables = doc.get("curves") or {}, doc.get("check_tables") or {}, {}
    files = {t["file"] for c in curves.values() for t in c["tables"]} | {c["file"] for c in checks.values()}
    for f in files:
        if not (base / f).exists():
            raise RegistryError(f"growth table {f} listed in curves.yaml does not exist")
        tables[f] = read_table(base / f)
    return GrowthData(curves, tables, checks)


def load_growth(data_dir: str | Path | None = None) -> GrowthData:
    return _load(str(Path(data_dir) if data_dir else DATA_DIR))
