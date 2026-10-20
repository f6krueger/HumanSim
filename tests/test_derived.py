"""Derived formulas. Several tests reproduce derived columns from the original parameter
document (placeholder numbers, but internally consistent), which checks the formulas."""

import pytest

from humansim import derived

# Projected-adult column of the original document (placeholder species, used only for formulas).
DOC = {
    "height": 215, "body_mass": 198.0, "arm_span": 230, "standing_reach": 280,
    "leg_length_trochanteric": 122, "arm_length": 97.5, "torso_length": 62,
    "upper_arm_length": 40.5, "forearm_length": 33.5, "femur_length": 53, "tibia_length": 54,
    "waist_circumference": 118, "hip_circumference": 140, "chest_circumference": 143,
    "head_height": 30, "penis_girth_erect": 15.0, "body_density": 1.19,
}


@pytest.fixture(scope="module")
def doc():
    return derived.compute(DOC)


def test_formulas_reproduce_document_columns(doc):
    assert doc["body_mass_index"] == pytest.approx(42.8, abs=0.05)
    assert doc["ape_index"] == pytest.approx(1.07, abs=0.005)
    assert doc["standing_reach_to_stature"] == pytest.approx(1.302, abs=0.001)
    assert doc["waist_to_hip_ratio"] == pytest.approx(0.84, abs=0.005)
    assert doc["chest_to_waist_ratio"] == pytest.approx(1.21, abs=0.005)
    assert doc["intermembral_index"] == pytest.approx(69.2, abs=0.05)
    assert doc["humerofemoral_index"] == pytest.approx(0.764, abs=0.001)
    assert doc["brachial_index"] == pytest.approx(0.827, abs=0.001)
    assert doc["crural_index"] == pytest.approx(1.019, abs=0.001)
    assert doc["leg_to_trunk_ratio"] == pytest.approx(1.97, abs=0.005)
    assert doc["arm_span_to_leg_length"] == pytest.approx(1.89, abs=0.005)
    assert doc["arm_length_to_leg_length"] == pytest.approx(0.799, abs=0.001)
    assert doc["penis_diameter_erect"] == pytest.approx(4.8, abs=0.05)
    assert doc["head_height_fraction"] == pytest.approx(0.140, abs=0.001)
    assert doc["body_volume"] == pytest.approx(166, abs=0.5)  # 198 kg / 1.19 kg/L


def test_bsa_reference_male():
    out = derived.compute({"body_mass": 73.0, "height": 176.0})
    assert out["body_surface_area"] == pytest.approx(1.89, abs=0.02)
    assert out["surface_to_mass"] == pytest.approx(out["body_surface_area"] / 73.0)


def test_physiology_formulas():
    out = derived.compute({"hemoglobin": 15.0, "sao2": 0.98, "pao2": 95.0,
                           "blood_volume": 5.3, "hematocrit": 0.45,
                           "cardiac_output": 6.5, "hepatic_flow_fraction": 0.255})
    assert out["arterial_o2_content"] == pytest.approx(200.0, abs=1.0)
    assert out["red_cell_volume"] == pytest.approx(2.385)
    assert out["hepatic_blood_flow"] == pytest.approx(1.6575)


def test_missing_inputs_give_none_not_zero():
    out = derived.compute({"body_mass": 73.0, "height": None})
    assert out["body_mass_index"] is None
    assert out["body_surface_area"] is None


def test_division_by_zero_gives_none():
    out = derived.compute({"body_mass": 73.0, "body_density": 0.0})
    assert out["body_volume"] is None


def test_every_derived_input_exists(registry):
    ids = set(registry.ids())
    for spec in derived.SPECS:
        assert set(spec.inputs) <= ids, spec.id
