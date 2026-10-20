import csv
from pathlib import Path

import pytest

from humansim import Context, load_couplings, load_modifiers, resolve, solve
from humansim.checks import bmi_consistency, growth_continuity
from humansim.growth import Z_P3, Z_P97, ADULT_AGE_YEARS, lms_value, load_growth
from humansim.registry import DATA_DIR

G = load_growth()


@pytest.fixture(scope="module")
def mods():
    return load_modifiers()


def rows(name):
    with open(Path(DATA_DIR) / "growth" / name, newline="") as fh:
        return list(csv.DictReader(fh))


# --- the tables and the LMS formula ----------------------------------------------------------------
@pytest.mark.parametrize("name", ["cdc2000_lenage.csv", "cdc2000_wtage.csv", "cdc2000_lenageinf.csv",
                                  "cdc2000_wtageinf.csv", "cdc2000_hcageinf.csv"])
def test_lms_formula_reproduces_the_published_percentiles(name):
    # Integrity check of both the files and the formula: P3, P50 and P97 must come back out of L, M, S.
    for r in rows(name)[::7]:
        L, M, S = float(r["L"]), float(r["M"]), float(r["S"])
        assert lms_value(L, M, S, 0) == pytest.approx(float(r["P50"]), rel=1e-6)
        assert lms_value(L, M, S, Z_P3) == pytest.approx(float(r["P3"]), rel=2e-3)
        assert lms_value(L, M, S, Z_P97) == pytest.approx(float(r["P97"]), rel=2e-3)


def test_known_cdc_values():
    assert G.at("height", "male", 0).value == pytest.approx(49.99, abs=0.01)
    assert G.at("body_mass", "male", 0).value == pytest.approx(3.53, abs=0.01)
    assert G.at("height", "female", 240).value == pytest.approx(163.34, abs=0.01)
    assert G.at("body_mass", "female", 240).value == pytest.approx(58.2, abs=0.05)


def test_interpolation_is_between_rows():
    a, b = G.at("height", "male", 60).value, G.at("height", "male", 60.5).value
    mid = G.at("height", "male", 60.25).value
    assert a < mid < b
    assert mid == pytest.approx(a + (b - a) * 0.5, rel=1e-3)  # 60.25 is halfway between rows 60 and 60.5


def test_out_of_table_age_gives_nothing():
    assert G.at("head_circumference", "male", 60) is None  # table ends at 36 months
    assert G.at("height", "male", 241) is None
    assert G.at("liver_mass", "male", 60) is None


def test_infant_to_child_table_jump_is_small():
    # recumbent length before 24 months, standing stature after: CDC splits the tables here
    a = G.at("height", "male", 23.99).value
    b = G.at("height", "male", 24.0).value
    assert abs(a - b) < 1.5


# --- pediatric context ----------------------------------------------------------------------------------
def test_child_gets_growth_values_only(registry, mods):
    s = resolve(registry, Context(sex="female", age=10), mods)
    assert s.get("height") == pytest.approx(G.at("height", "female", 120).value)
    assert {"growth_curve", "interpolated"} <= s.flags["height"]
    # no adult value leaks into a child
    for pid in ("arm_span", "hematocrit", "hand_length", "neck_circumference"):
        assert s.get(pid) is None, pid
        assert "missing" in s.flags[pid]
    assert s.get("body_mass_index") == pytest.approx(s.get("body_mass") / (s.get("height") / 100) ** 2)
    assert "interpolated" in s.flags["body_mass_index"]
    assert any("Pediatric" in n for n in s.notes)


def test_sex_difference_emerges_from_the_curves(registry, mods):
    m = resolve(registry, Context(sex="male", age=18), mods)
    f = resolve(registry, Context(sex="female", age=18), mods)
    assert m.get("height") > f.get("height") and m.get("body_mass") > f.get("body_mass")


def test_growth_z_moves_the_value_and_the_range_stays(registry, mods):
    med = resolve(registry, Context(age=8), mods)
    tall = resolve(registry, Context(age=8, growth_z=1.0), mods)
    assert tall.get("height") > med.get("height")
    assert tall.ranges["height"] == med.ranges["height"]
    assert med.ranges["height"][0] < med.get("height") < med.ranges["height"][1]


def test_pediatric_ranges_replace_adult_ranges(registry, mods):
    low, high = resolve(registry, Context(age=8), mods).ranges["height"]
    # 150 cm is inside the adult female range but far above the 97th percentile of an 8-year-old
    s = resolve(registry, Context(sex="female", age=8, edits={"height": 150.0}), mods)
    assert 150 > high and "outside_reference_range" in s.flags["height"]
    # an edited adult-only parameter in a child is not judged against the adult range
    s = resolve(registry, Context(age=8, edits={"neck_circumference": 5.0}), mods)
    assert "outside_reference_range" not in s.flags["neck_circumference"]


def test_age_20_and_over_is_the_adult_baseline_with_a_note(registry, mods):
    adult = resolve(registry, Context(), mods)
    for age in (20, 35, 60):
        s = resolve(registry, Context(age=age), mods)
        assert s.get("height") == adult.get("height") and s.ranges is None
    assert not resolve(registry, Context(age=20), mods).notes
    assert any("ageing" in n for n in resolve(registry, Context(age=60), mods).notes)
    assert resolve(registry, Context(age=19.9), mods).ranges is not None


def test_bad_age_is_rejected(registry, mods):
    with pytest.raises(ValueError, match="age"):
        resolve(registry, Context(age=-1), mods)


# --- checks ------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("sex", ["male", "female"])
@pytest.mark.parametrize("age", [2, 5, 10, 15, 19])
def test_bmi_from_medians_matches_cdc_bmi_median(registry, mods, sex, age):
    r = bmi_consistency(resolve(registry, Context(sex=sex, age=age), mods))
    assert r.status == "pass", r.message


def test_continuity_with_adult_is_reported_honestly(registry):
    res = {r.id: r for r in growth_continuity(registry, "male")}
    assert res["growth_continuity:height:male"].status == "pass"
    # the Army sample is much heavier than the general population at 20: a warning, not a hidden gap
    assert res["growth_continuity:body_mass:male"].status == "warn"


# --- solver at a pediatric age -------------------------------------------------------------------------------------
def test_solver_works_with_what_a_child_has(registry, mods):
    cps = load_couplings(registry)
    s = solve(registry, Context(age=3, edits={"height": 100.0}), mods, cps)
    assert s.solver["converged"] and s.get("height") == 100.0
    base = resolve(registry, Context(age=3), mods)
    assert s.get("body_mass") > base.get("body_mass")  # taller child is heavier
    assert s.get("liver_mass") == pytest.approx(base.get("liver_mass"), rel=0.3)
    assert s.ranges is not None


def test_unedited_child_is_a_solver_fixed_point(registry, mods):
    s = solve(registry, Context(sex="female", age=12), mods, load_couplings(registry))
    assert s.changed_from(resolve(registry, Context(sex="female", age=12), mods)) == {}


# --- ICRP anchor ages -------------------------------------------------------------------------------------
def test_icrp_anchor_age_is_exact_and_not_interpolated(registry, mods):
    s = resolve(registry, Context(sex="female", age=10), mods)
    assert s.get("liver_mass") == 0.83 and "interpolated" not in s.flags["liver_mass"]
    assert s.get("cardiac_output") == 5.0
    assert resolve(registry, Context(age=0), mods).get("liver_mass") == 0.13


def test_between_anchors_is_linear_and_flagged(registry, mods):
    s = resolve(registry, Context(age=7.5), mods)  # halfway between 5 y (0.57 kg) and 10 y (0.83 kg)
    assert s.get("liver_mass") == pytest.approx(0.70)
    assert "interpolated" in s.flags["liver_mass"]
    assert "interpolated" in s.flags["hepatic_blood_flow"] or s.get("hepatic_blood_flow") is None
    # 15 y male to adult: interpolated towards the adult value at age 20
    m = resolve(registry, Context(age=17.5), mods)
    assert m.get("liver_mass") == pytest.approx((1.3 + 1.8) / 2)


def test_adult_only_values_do_not_reach_a_child(registry, mods):
    s = resolve(registry, Context(age=10), mods)
    assert s.get("hepatic_flow_fraction") is None  # adult-only in the report
    assert s.get("hematocrit") is None


def test_sex_specific_organs(registry, mods):
    assert resolve(registry, Context(sex="female", age="adult"), mods).get("testes_mass") is None
    assert resolve(registry, Context(sex="male", age="adult"), mods).get("ovaries_mass") is None
    assert resolve(registry, Context(sex="male", age="adult"), mods).get("testes_mass") == 0.035
    assert resolve(registry, Context(sex="female", age=15), mods).get("uterus_mass") == 0.03
