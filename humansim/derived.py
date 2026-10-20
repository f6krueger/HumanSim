"""Derived layer. Derived values are never stored; they are computed from anchors.

Each spec lists its inputs by parameter id. A derived value is computed only when every
input has a value. Derived values may depend on other derived values.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Mapping


@dataclass(frozen=True)
class DerivedSpec:
    id: str
    name: str
    unit: str
    group: str
    inputs: tuple[str, ...]
    fn: Callable[..., float]
    formula: str
    source_ref: str | None = None
    note: str = ""


SPECS: list[DerivedSpec] = []


def _add(id, name, unit, group, inputs, fn, formula, source_ref=None, note=""):
    SPECS.append(DerivedSpec(id, name, unit, group, tuple(inputs), fn, formula, source_ref, note))


def _ratio(num: str, den: str) -> Callable[..., float]:
    return lambda **kw: kw[num] / kw[den]


# --- body size and surface -------------------------------------------------------------
_add("body_mass_index", "Body mass index", "kg/m^2", "anthropometry",
     ("body_mass", "height"), lambda body_mass, height: body_mass / (height / 100) ** 2,
     "mass_kg / height_m^2")
_add("body_surface_area", "Body surface area", "m^2", "anthropometry",
     ("body_mass", "height"),
     lambda body_mass, height: 0.007184 * body_mass ** 0.425 * height ** 0.725,
     "0.007184 * mass_kg^0.425 * height_cm^0.725", source_ref="dubois1916",
     note="Du Bois formula. Other formulas (Mosteller etc.) give slightly different values.")
_add("surface_to_mass", "Surface-to-mass ratio", "m^2/kg", "anthropometry",
     ("body_surface_area", "body_mass"), _ratio("body_surface_area", "body_mass"),
     "body_surface_area / body_mass")
_add("body_volume", "Body volume", "L", "anthropometry",
     ("body_mass", "body_density"), _ratio("body_mass", "body_density"),
     "body_mass_kg / density_(g/cm^3)", note="1 g/cm^3 equals 1 kg/L.")

# --- proportions -----------------------------------------------------------------------
_add("ape_index", "Ape index (arm span / stature)", "dimensionless", "proportions",
     ("arm_span", "height"), _ratio("arm_span", "height"), "arm_span / height")
_add("waist_to_hip_ratio", "Waist-to-hip ratio", "dimensionless", "proportions",
     ("waist_circumference", "hip_circumference"),
     _ratio("waist_circumference", "hip_circumference"), "waist / hip")
_add("chest_to_waist_ratio", "Chest-to-waist ratio", "dimensionless", "proportions",
     ("chest_circumference", "waist_circumference"),
     _ratio("chest_circumference", "waist_circumference"), "chest / waist")
_add("intermembral_index", "Intermembral index", "dimensionless", "proportions",
     ("upper_arm_length", "forearm_length", "femur_length", "tibia_length"),
     lambda upper_arm_length, forearm_length, femur_length, tibia_length:
     100 * (upper_arm_length + forearm_length) / (femur_length + tibia_length),
     "100 * (upper_arm + forearm) / (femur + tibia)")
_add("humerofemoral_index", "Humerofemoral index", "dimensionless", "proportions",
     ("upper_arm_length", "femur_length"), _ratio("upper_arm_length", "femur_length"),
     "upper_arm / femur")
_add("brachial_index", "Brachial index", "dimensionless", "proportions",
     ("forearm_length", "upper_arm_length"), _ratio("forearm_length", "upper_arm_length"),
     "forearm / upper_arm")
_add("crural_index", "Crural index", "dimensionless", "proportions",
     ("tibia_length", "femur_length"), _ratio("tibia_length", "femur_length"),
     "tibia / femur")
_add("leg_to_trunk_ratio", "Leg length / trunk length", "dimensionless", "proportions",
     ("leg_length_trochanteric", "torso_length"),
     _ratio("leg_length_trochanteric", "torso_length"), "trochanteric_leg / torso")
_add("arm_span_to_leg_length", "Arm span / leg length", "dimensionless", "proportions",
     ("arm_span", "leg_length_trochanteric"), _ratio("arm_span", "leg_length_trochanteric"),
     "arm_span / trochanteric_leg")
_add("arm_length_to_leg_length", "Arm length / leg length", "dimensionless", "proportions",
     ("arm_length", "leg_length_trochanteric"), _ratio("arm_length", "leg_length_trochanteric"),
     "arm_length / trochanteric_leg")
_add("standing_reach_to_stature", "Standing reach / stature", "dimensionless", "proportions",
     ("standing_reach", "height"), _ratio("standing_reach", "height"),
     "standing_reach / height")
_add("penis_diameter_erect", "Penis diameter, erect", "cm", "anthropometry",
     ("penis_girth_erect",), lambda penis_girth_erect: penis_girth_erect / math.pi,
     "girth / pi", note="Assumes a circular cross-section.")

SEGMENTS = {
    "head_height": "Head height", "head_length": "Head length", "head_breadth": "Head breadth",
    "neck_length": "Neck length", "shoulder_height": "Shoulder height",
    "torso_length": "Torso length", "femur_length": "Femur length",
    "tibia_length": "Tibia length", "foot_length": "Foot length",
    "foot_breadth": "Foot breadth", "ankle_height": "Ankle height",
    "upper_arm_length": "Upper arm length", "forearm_length": "Forearm length",
    "hand_length": "Hand length", "hand_breadth": "Hand breadth",
}
for _sid, _sname in SEGMENTS.items():
    _add(f"{_sid}_fraction", f"{_sname}, fraction of stature", "dimensionless", "proportions",
         (_sid, "height"), _ratio(_sid, "height"), f"{_sid} / height")

# --- composition -----------------------------------------------------------------------
_add("muscle_fraction", "Skeletal muscle fraction of body mass", "dimensionless", "composition",
     ("skeletal_muscle_mass", "body_mass"), _ratio("skeletal_muscle_mass", "body_mass"),
     "skeletal_muscle_mass / body_mass")
_add("skeleton_fraction", "Skeleton fraction of body mass", "dimensionless", "composition",
     ("skeleton_mass", "body_mass"), _ratio("skeleton_mass", "body_mass"),
     "skeleton_mass / body_mass")
_add("adipose_fraction", "Adipose tissue fraction of body mass", "dimensionless", "composition",
     ("adipose_tissue_mass", "body_mass"), _ratio("adipose_tissue_mass", "body_mass"),
     "adipose_tissue_mass / body_mass",
     note="Adipose tissue is not the same as body fat: adipose tissue includes water and protein.")

# --- cardiovascular and haematological --------------------------------------------------
_add("red_cell_volume", "Red cell volume", "L", "physiology",
     ("blood_volume", "hematocrit"),
     lambda blood_volume, hematocrit: blood_volume * hematocrit,
     "blood_volume * hematocrit")
_add("arterial_o2_content", "Arterial oxygen content", "mL/L", "physiology",
     ("hemoglobin", "sao2", "pao2"),
     lambda hemoglobin, sao2, pao2: 10 * (1.34 * hemoglobin * sao2 + 0.0031 * pao2),
     "10 * (1.34 * Hb_g/dL * SaO2 + 0.0031 * PaO2_mmHg)",
     note="Standard clinical constants (1.34 mL O2 per g Hb; 0.0031 mL/dL/mmHg dissolved). "
          "Recalled, not retrieved.")
_add("hepatic_blood_flow", "Hepatic blood flow (total)", "L/min", "physiology",
     ("cardiac_output", "hepatic_flow_fraction"),
     lambda cardiac_output, hepatic_flow_fraction: cardiac_output * hepatic_flow_fraction,
     "cardiac_output * hepatic_flow_fraction")


def compute(values: Mapping[str, float | None]) -> dict[str, float | None]:
    """Compute every derived value that can be computed from `values`."""
    known: dict[str, float | None] = dict(values)
    out: dict[str, float | None] = {}
    pending = list(SPECS)
    progress = True
    while progress and pending:
        progress = False
        for spec in list(pending):
            if all(known.get(i) is not None for i in spec.inputs):
                try:
                    v: float | None = float(spec.fn(**{i: known[i] for i in spec.inputs}))
                    if math.isnan(v) or math.isinf(v):
                        v = None
                except (ZeroDivisionError, ValueError, OverflowError):
                    v = None
                known[spec.id] = out[spec.id] = v
                pending.remove(spec)
                progress = True
    for spec in pending:
        out[spec.id] = None
    return out
