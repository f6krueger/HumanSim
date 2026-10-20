"""Unit handling. All unit strings in the data files are validated with pint."""

from __future__ import annotations

import pint

ureg = pint.UnitRegistry()
# Biological activity has no fixed mass equivalent, so it gets its own dimension (not convertible to mass).
ureg.define("international_unit = [biological_activity] = IU")


class UnitError(ValueError):
    """Raised for an unknown unit or a value that cannot be converted."""


def validate_unit(unit: str) -> None:
    try:
        ureg.Unit(unit)
    except Exception as exc:  # pint raises several exception types
        raise UnitError(f"Unknown unit {unit!r}: {exc}") from exc


def convert(value: float, from_unit: str, to_unit: str) -> float:
    try:
        return float(ureg.Quantity(value, from_unit).to(to_unit).magnitude)
    except Exception as exc:
        raise UnitError(f"Cannot convert {value} {from_unit} to {to_unit}: {exc}") from exc


def parse_value(text: str, target_unit: str) -> float:
    """Parse '90' (taken to be in target_unit) or '198 lb' (converted to target_unit)."""
    text = text.strip()
    try:
        return float(text)
    except ValueError:
        pass
    try:
        return float(ureg.Quantity(text).to(target_unit).magnitude)
    except Exception as exc:
        raise UnitError(f"Cannot read {text!r} as a value in {target_unit}: {exc}") from exc
