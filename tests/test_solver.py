import pytest

from humansim import Context, Registry, load_couplings, load_modifiers, resolve, solve
from humansim.checks import coupling_checks, rmr_organ_crosscheck
from humansim.couplings import Coupling
from humansim.registry import RegistryError


@pytest.fixture(scope="module")
def mods():
    return load_modifiers()


@pytest.fixture(scope="module")
def couplings(registry):
    return load_couplings(registry)


def run(registry, mods, couplings, sex="male", **edits):
    return solve(registry, Context(sex=sex, edits=edits), mods, couplings)


def changed(registry, mods, state, sex="male"):
    return set(state.changed_from(resolve(registry, Context(sex=sex), mods)))


# --- the unedited model is a fixed point ------------------------------------------------------
@pytest.mark.parametrize("sex", ["male", "female"])
def test_no_edits_changes_nothing(registry, mods, couplings, sex):
    s = run(registry, mods, couplings, sex)
    assert s.solver["converged"]
    assert changed(registry, mods, s, sex) == set()
    assert s.trace == [t for t in resolve(registry, Context(sex=sex), mods).trace]


# --- deadband: normal variation moves nothing, essential change adjusts ------------------------
def test_small_liver_change_is_inside_the_deadband(registry, mods, couplings):
    s = run(registry, mods, couplings, liver_mass=1.9)  # +5.6%
    moved = changed(registry, mods, s)
    assert "cardiac_output" not in moved and "heart_mass" not in moved
    assert "resting_metabolic_rate" in moved  # accounting: no deadband
    assert s.get("cardiac_output") == 6.5


def test_mass_within_reference_range_leaves_organs_alone_but_outside_scales_them(registry, mods, couplings):
    inside = run(registry, mods, couplings, body_mass=100.0)  # inside 60.8-116.4
    assert inside.get("liver_mass") == 1.8 and inside.get("skeletal_muscle_mass") == 29.0
    outside = run(registry, mods, couplings, body_mass=130.0)
    assert outside.get("liver_mass") > 1.8
    assert "solved" in outside.flags["liver_mass"] and "coupling_unsourced" in outside.flags["liver_mass"]


def test_response_is_continuous_at_the_band_edge(registry, mods, couplings):
    hi = registry["body_mass"].ref_range["male"][1]
    just_out = run(registry, mods, couplings, body_mass=hi * 1.0001)
    assert just_out.get("liver_mass") == pytest.approx(1.8, rel=1e-3)


# --- the liver chain ----------------------------------------------------------------------------
def test_liver_nine_kg_chain(registry, mods, couplings):
    s = run(registry, mods, couplings, liver_mass=9.0)
    assert s.solver["converged"]
    assert s.get("liver_mass") == 9.0 and s.flags["liver_mass"] == {"edited"}
    assert s.get("body_mass") == pytest.approx(85.52 + 7.2)  # the extra liver is extra mass
    assert s.get("resting_metabolic_rate") > 1.7 * resolve(registry, Context(), mods).get("resting_metabolic_rate")
    assert s.get("cardiac_output") > 6.5 and s.get("heart_mass") > 0.33
    assert s.get("hepatic_flow_fraction") > 0.255
    assert s.get("hepatic_blood_flow") == pytest.approx(s.get("cardiac_output") * s.get("hepatic_flow_fraction"))
    # five times the reference liver is far outside what the relationships were built for
    assert "outside_coupling_validity" in s.flags["hepatic_flow_fraction"]
    assert [t.source for t in s.trace if t.stage == "edit"] == ["user"]


def test_rmr_matches_organ_sum_when_organs_change(registry, mods, couplings):
    s = run(registry, mods, couplings, liver_mass=9.0)
    assert abs(rmr_organ_crosscheck(s).details["difference"]) < 0.20  # body_mass basis after an edit


def test_edited_target_is_never_overwritten(registry, mods, couplings):
    s = run(registry, mods, couplings, liver_mass=9.0, cardiac_output=6.5)
    assert s.get("cardiac_output") == 6.5 and s.flags["cardiac_output"] == {"edited"}
    assert "cardiac_output_follows_oxygen_demand" in s.solver["skipped"]


# --- Fick principle -------------------------------------------------------------------------------
def test_anaemia_raises_cardiac_output_and_heart_mass(registry, mods, couplings):
    s = run(registry, mods, couplings, "female", hemoglobin=9.0)
    assert s.get("arterial_o2_content") < 130
    assert s.get("cardiac_output") > 5.9 * 1.3
    assert s.get("heart_mass") > 0.25


def test_mild_haemoglobin_change_is_inside_the_deadband(registry, mods, couplings):
    assert run(registry, mods, couplings, hemoglobin=14.5).get("cardiac_output") == 6.5


# --- height drives shape and mass -----------------------------------------------------------------
def test_tall_male_scales_dimensions_mass_and_head(registry, mods, couplings):
    s = run(registry, mods, couplings, height=213.36)  # 7 ft
    k = 213.36 / 175.62
    assert s.get("arm_span") == pytest.approx(181.42 * k, rel=1e-3)
    assert s.get("body_mass") == pytest.approx(85.52 * k ** 2, rel=1e-3)
    assert s.get("head_circumference") == pytest.approx(57.44 * k ** 0.5, rel=1e-3)
    assert s.get("head_circumference") / 57.44 < s.get("arm_span") / 181.42  # head grows less
    assert s.get("body_mass_index") == pytest.approx(resolve(registry, Context(), mods).get("body_mass_index"), rel=1e-3)
    assert "solved" in s.flags["arm_span"] and "outside_reference_range" in s.flags["height"]


# --- female goes through the same couplings ---------------------------------------------------------
def test_female_chain_uses_female_reference(registry, mods, couplings):
    s = run(registry, mods, couplings, "female", liver_mass=2.8)
    assert s.get("body_mass") == pytest.approx(67.76 + 1.4)
    assert s.get("cardiac_output") == pytest.approx(5.9, rel=0.01) or s.get("cardiac_output") > 5.9
    assert s.get("hepatic_flow_fraction") > 0.255


def test_missing_reference_skips_coupling_and_says_so(registry, mods, couplings):
    s = run(registry, mods, couplings, height=200.0)
    assert s.get("femur_length") is None
    assert "linear_dimensions_follow_height:femur_length" in s.solver["skipped"]


# --- loader and checks ---------------------------------------------------------------------------------
def test_select_and_targets_expand(registry, couplings):
    ids = {c.id for c in couplings}
    assert "linear_dimensions_follow_height:arm_span" in ids
    assert "linear_dimensions_follow_height:height" not in ids
    assert "head_dimensions_follow_height:head_length" in ids
    assert coupling_checks(registry, couplings)[0].status == "pass"


def test_bad_couplings_are_rejected_or_reported(registry):
    with pytest.raises(RegistryError):
        Coupling("a", "x", "nope")
    with pytest.raises(RegistryError):
        Coupling("a", "x", "power_ratio", role="accounting", deadband=0.1)
    a = Coupling("a", "body_mass_index", "power_ratio", args={"inputs": {"height": 1}})
    b = Coupling("b", "liver_mass", "power_ratio", args={"inputs": {"nope": 1}})
    msg = coupling_checks(registry, [a, b])[0]
    assert msg.status == "fail" and "derived" in msg.message and "unknown input" in msg.message


def test_feedback_loop_that_diverges_is_reported(registry, mods):
    # heart mass follows cardiac output with exponent 6: loop gain above 1, so no fixed point
    cps = [Coupling("rmr", "resting_metabolic_rate", "rmr_organ_sum", kind="mechanistic", role="accounting",
                    args={"total": "body_mass", "residual_rate": 12.0,
                          "rates": {"heart_mass": 4400.0, "skeletal_muscle_mass": 13.0}}),
           Coupling("co", "cardiac_output", "fick_cardiac_output", args={"demand": "resting_metabolic_rate"}),
           Coupling("h", "heart_mass", "power_ratio", args={"inputs": {"cardiac_output": 6.0}})]
    s = solve(registry, Context(edits={"heart_mass": 0.33 * 1.0}), mods, cps)  # heart pinned: no loop
    assert s.solver["converged"]
    s = solve(registry, Context(edits={"skeletal_muscle_mass": 60.0}), mods, cps, max_iter=50)
    assert s.solver["converged"] is False
