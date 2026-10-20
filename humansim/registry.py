"""Parameter registry: loads value files, validates them, and serves lookups.

Anchors, regulated and dynamic parameters live in YAML files under data/values/.
Derived parameters are defined in code (derived.py) because they carry formulas.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Mapping

import yaml

from .units import validate_unit

DATA_DIR = Path(__file__).parent / "data"

KINDS = ("anchor", "derived", "regulated", "dynamic")
# "missing" means the parameter is defined but has no value yet.
STATUSES = ("referenced", "unsourced", "placeholder", "missing")

_ID_RE = re.compile(r"^[a-z][a-z0-9_]*$")
_ALLOWED_KEYS = {
    "id", "name", "unit", "kind", "values", "sd", "status", "source_ref", "note",
    "ref_range", "range_source_ref", "viability_limit",
}


class RegistryError(ValueError):
    """Raised when a data file is malformed."""


@dataclass(frozen=True)
class Parameter:
    id: str
    name: str
    unit: str
    kind: str
    group: str
    values: Mapping[str, Mapping[object, float]] = field(default_factory=dict)
    sd: Mapping[str, Mapping[object, float]] = field(default_factory=dict)
    status: str = "missing"
    source_ref: str | None = None
    note: str = ""
    ref_range: tuple[float, float] | None = None
    range_source_ref: str | None = None
    viability_limit: tuple[float, float] | None = None
    formula: str | None = None

    def baseline(self, sex: str, age: object) -> float | None:
        """Mean value for (sex, age), or None if the data has no entry for it."""
        return self.values.get(sex, {}).get(age)

    def spread(self, sex: str, age: object) -> float | None:
        """Standard deviation for (sex, age) if the data has one."""
        return self.sd.get(sex, {}).get(age)


class Registry:
    def __init__(self, parameters: Mapping[str, Parameter]):
        self._p = dict(parameters)

    def __getitem__(self, pid: str) -> Parameter:
        try:
            return self._p[pid]
        except KeyError:
            raise KeyError(f"Unknown parameter {pid!r}") from None

    def __contains__(self, pid: object) -> bool:
        return pid in self._p

    def __iter__(self) -> Iterator[Parameter]:
        return iter(self._p.values())

    def __len__(self) -> int:
        return len(self._p)

    def ids(self) -> list[str]:
        return list(self._p)

    def by_kind(self, kind: str) -> list[Parameter]:
        return [p for p in self if p.kind == kind]

    def by_group(self, group: str) -> list[Parameter]:
        return [p for p in self if p.group == group]

    def groups(self) -> list[str]:
        seen: list[str] = []
        for p in self:
            if p.group not in seen:
                seen.append(p.group)
        return seen

    @classmethod
    def load(cls, data_dir: str | Path | None = None) -> "Registry":
        root = Path(data_dir) if data_dir else DATA_DIR
        files = sorted((root / "values").glob("*.yaml"))
        if not files:
            raise RegistryError(f"No value files found in {root / 'values'}")
        params: dict[str, Parameter] = {}
        for path in files:
            doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            if "group" not in doc or "parameters" not in doc:
                raise RegistryError(f"{path.name}: needs 'group' and 'parameters'")
            for raw in doc["parameters"]:
                p = _build(raw, doc["group"], path.name)
                if p.id in params:
                    raise RegistryError(f"{path.name}: duplicate parameter id {p.id!r}")
                params[p.id] = p

        from .derived import SPECS  # local import: derived.py has no registry dependency

        for spec in SPECS:
            if spec.id in params:
                raise RegistryError(f"Derived id {spec.id!r} collides with a stored parameter")
            validate_unit(spec.unit)
            for i in spec.inputs:
                if i not in params and i not in {s.id for s in SPECS}:
                    raise RegistryError(f"Derived {spec.id!r} uses unknown input {i!r}")
            params[spec.id] = Parameter(
                id=spec.id, name=spec.name, unit=spec.unit, kind="derived", group=spec.group,
                status="derived", source_ref=spec.source_ref, note=spec.note, formula=spec.formula,
            )
        return cls(params)


def _range(raw: object, where: str) -> tuple[float, float] | None:
    if raw is None:
        return None
    if not (isinstance(raw, (list, tuple)) and len(raw) == 2):
        raise RegistryError(f"{where}: range must be [low, high]")
    lo, hi = float(raw[0]), float(raw[1])
    if lo >= hi:
        raise RegistryError(f"{where}: range low must be below high")
    return (lo, hi)


def _build(raw: dict, group: str, filename: str) -> Parameter:
    pid = raw.get("id", "?")
    where = f"{filename}:{pid}"
    unknown = set(raw) - _ALLOWED_KEYS
    if unknown:
        raise RegistryError(f"{where}: unknown keys {sorted(unknown)}")
    for key in ("id", "name", "unit", "kind"):
        if key not in raw:
            raise RegistryError(f"{where}: missing required key {key!r}")
    if not _ID_RE.match(pid):
        raise RegistryError(f"{where}: id must be snake_case")
    if raw["kind"] not in KINDS or raw["kind"] == "derived":
        raise RegistryError(f"{where}: kind must be anchor, regulated or dynamic")
    try:
        validate_unit(raw["unit"])
    except ValueError as exc:
        raise RegistryError(f"{where}: {exc}") from exc

    values: dict[str, dict[object, float]] = {}
    for sex, by_age in (raw.get("values") or {}).items():
        if sex not in ("male", "female"):
            raise RegistryError(f"{where}: sex must be male or female, got {sex!r}")
        values[sex] = {age: float(v) for age, v in by_age.items()}

    sd: dict[str, dict[object, float]] = {}
    for sex, by_age in (raw.get("sd") or {}).items():
        sd[sex] = {}
        for age, s in by_age.items():
            if float(s) <= 0:
                raise RegistryError(f"{where}: sd must be positive")
            if age not in values.get(sex, {}):
                raise RegistryError(f"{where}: sd given for ({sex}, {age}) but there is no value")
            sd[sex][age] = float(s)

    status = raw.get("status", "missing" if not values else "unsourced")
    if status not in STATUSES:
        raise RegistryError(f"{where}: status must be one of {STATUSES}")
    if (status == "missing") != (not values):
        raise RegistryError(f"{where}: status is 'missing' exactly when there are no values")

    return Parameter(
        id=pid, name=raw["name"], unit=raw["unit"], kind=raw["kind"], group=group,
        values=values, sd=sd, status=status, source_ref=raw.get("source_ref"),
        note=raw.get("note", ""), ref_range=_range(raw.get("ref_range"), where),
        range_source_ref=raw.get("range_source_ref"),
        viability_limit=_range(raw.get("viability_limit"), where),
    )
