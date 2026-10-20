import pytest

from humansim import Context, Modifier, Parameter, Registry, load_modifiers, resolve
from humansim.checks import female_rule_coverage, modifier_checks, range_flags
from humansim.validation import (classify_ratio, load_validation, scaling_analysis,
                                 scaling_expectation)


@pytest.fixture(scope="module")
def mods():
    return load_modifiers()


@pytest.fixture(scope="module")
def female_doc():
    return load_validation("ansur2_female")


@pytest.fixture(scope="module")
def female(registry, mods):
    return resolve(registry, Context(sex="female"), mods)


@pytest.fixture(scope="module")
def male(registry, mods):
    return resolve(registry, Context(sex="male"), mods)


# --- calibration (consistency, not independent validation) --------------------------------
def test_derived_female_reproduces_ansur_female_means(female, female_doc):
    for pid, rec in female_doc["parameters"].items():
        assert female.get(pid) == pytest.approx(rec["mean"], rel=1e-3), pid


def test_male_is_unchanged_by_sex_rules(registry, mods):
    with_rules = resolve(registry, Context(sex="male"), mods)
    without = resolve(registry, Context(sex="male"), [])
    assert with_rules.values == without.values
    assert with_rules.trace == []


def test_every_sex_rule_targets_a_valid_parameter_once(registry, mods):
    assert modifier_checks(registry, mods)[0].status == "pass"
    assert len(mods) == 38
    assert {m.stage for m in mods} == {"sex"}
    assert all(m.when == {"sex": "female"} for m in mods)


def test_rule_types_are_labelled_and_sourced(mods):
    types = {m.rule_type for m in mods}
    assert types == {"scaling", "anatomical", "physiological", "assumed_neutral"}
    for m in mods:
        assert m.status in ("referenced", "unsourced")
        if m.status == "referenced":
            assert m.source_ref
        if m.rule_type == "assumed_neutral":
            assert m.status == "unsourced"


# --- the ICRP female/male pairs -------------------------------------------------------------
def test_icrp_female_values(female):
    assert female.get("hematocrit") == pytest.approx(0.42, abs=1e-3)
    assert female.get("cardiac_output") == pytest.approx(5.9, abs=1e-3)
    assert female.get("blood_volume") == pytest.approx(3.9, abs=1e-3)


# --- parameters without a rule are not silently filled with male values ---------------------
def test_no_rule_means_no_value_for_female(female, male):
    assert male.get("liver_mass") is not None
    assert female.get("liver_mass") is None
    assert female.flags["liver_mass"] == {"missing", "no_sex_rule"}
    assert female.get("hepatic_blood_flow") is None  # derived from a parameter with no female value


def test_keep_unadjusted_keeps_the_male_value_but_flags_it(registry, mods, male):
    s = resolve(registry, Context(sex="female", keep_unadjusted=True), mods)
    assert s.get("liver_mass") == male.get("liver_mass")
    assert {"derived_from_male", "no_sex_rule"} <= s.flags["liver_mass"]
    assert "no_sex_rule" in s.flags["hepatic_blood_flow"] or s.get("hepatic_blood_flow") is not None


def test_assumed_neutral_parameters_are_kept_and_flagged_unsourced(female):
    assert female.get("core_temperature") == 37.0
    assert "modifier_unsourced" in female.flags["core_temperature"]
    assert "no_sex_rule" not in female.flags["core_temperature"]


def test_female_derived_values_follow_from_derived_inputs(female, male):
    h, m = female.get("height"), female.get("body_mass")
    assert female.get("body_mass_index") == pytest.approx(m / (h / 100) ** 2)
    assert female.get("body_mass_index") < male.get("body_mass_index")
    assert female.get("waist_to_hip_ratio") < male.get("waist_to_hip_ratio")


def test_female_edits_override_and_propagate(registry, mods):
    s = resolve(registry, Context(sex="female", edits={"body_mass": 80.0, "liver_mass": 1.5}), mods)
    assert s.get("body_mass") == 80.0
    assert s.get("liver_mass") == 1.5 and s.flags["liver_mass"] == {"edited"}
    assert [t.stage for t in s.trace if t.param == "body_mass"] == ["sex", "edit"]


def test_female_coverage_check_lists_missing_rules(female):
    r = female_rule_coverage(female)
    assert r.status == "warn" and "liver_mass" in r.details["no_sex_rule"]
    assert female_rule_coverage(resolve(Registry.load(), Context())).status == "skipped"


# --- reference ranges are per sex -------------------------------------------------------------
def test_female_baseline_is_inside_female_ranges_not_flagged(female):
    assert range_flags(female).status == "pass"
    assert "outside_reference_range" not in female.flags["neck_circumference"]


def test_male_range_is_not_applied_to_a_female(registry, mods):
    # 185 cm is inside the male range (up to 189.7) but far above the female 97.5th percentile (175.8).
    assert "outside_reference_range" not in resolve(registry, Context(edits={"height": 185.0}), mods).flags["height"]
    assert "outside_reference_range" in resolve(registry, Context(sex="female", edits={"height": 185.0}), mods).flags["height"]
    # Haemoglobin only has a male range, so a female value is never flagged against it.
    s = resolve(registry, Context(sex="female", edits={"hemoglobin": 25.0}), mods)
    assert "outside_reference_range" not in s.flags["hemoglobin"]


# --- stored female values take precedence over rules ------------------------------------------
def test_stored_female_value_is_not_modified_twice():
    a = Parameter(id="a", name="A", unit="kg", kind="anchor", group="g", status="referenced",
                  source_ref="s", values={"male": {"adult": 100.0}, "female": {"adult": 60.0}})
    reg = Registry({"a": a})
    rule = Modifier("r", "sex", "ratio", {"a": 0.5}, when={"sex": "female"})
    s = resolve(reg, Context(sex="female"), [rule])
    assert s.get("a") == 60.0
    assert "stored female value" in s.trace[0].note
    assert resolve(reg, Context(sex="male"), [rule]).get("a") == 100.0


def test_duplicate_sex_rules_are_reported(registry):
    a = Modifier("a", "sex", "ratio", {"height": 0.9}, when={"sex": "female"})
    b = Modifier("b", "sex", "ratio", {"height": 0.8}, when={"sex": "female"})
    r = modifier_checks(registry, [a, b])[0]
    assert r.status == "fail" and "duplicate" in r.message
    bad = modifier_checks(registry, [Modifier("c", "sex", "ratio", {"body_mass_index": 1.0})])[0]
    assert bad.status == "fail" and "derived" in bad.message


# --- scaling versus anatomical -------------------------------------------------------------------
def test_scaling_expectation_and_classification():
    assert scaling_expectation("cm", 0.93) == 0.93
    assert scaling_expectation("kg", 0.9) == pytest.approx(0.729)
    assert scaling_expectation("mL", 0.9) is None
    assert classify_ratio(0.94, 0.93) == "scaling"
    assert classify_ratio(0.85, 0.93) == "anatomical"
    assert classify_ratio(1.0, None) == "anatomical"


def test_known_dimorphisms_are_classified_anatomical(registry, female_doc):
    rows = {r["id"]: r for r in scaling_analysis(female_doc, registry)}
    # Narrower shoulders, thinner neck and arms, wider pelvis and hips than size alone predicts.
    for pid in ("biacromial_breadth", "neck_circumference", "upper_arm_circ_flexed"):
        assert rows[pid]["type"] == "anatomical" and rows[pid]["deviation"] < 0, pid
    for pid in ("bi_iliac_breadth", "hip_circumference", "thigh_circumference"):
        assert rows[pid]["type"] == "anatomical" and rows[pid]["deviation"] > 0, pid
    for pid in ("height", "sitting_height", "crotch_height"):
        assert rows[pid]["type"] == "scaling", pid


def test_cli_sex_commands():
    import subprocess, sys

    def run(*a):
        return subprocess.run([sys.executable, "-m", "humansim", *a], capture_output=True, text=True, timeout=60)

    out = run("compare")
    assert out.returncode == 0 and "F/M" in out.stdout and "hematocrit" in out.stdout
    out = run("validate")
    assert out.returncode == 0 and "NOT independent validation" in out.stdout
    out = run("show", "--sex", "female", "--group", "anthropometry", "--trace")
    assert out.returncode == 0 and "sex_female_height" in out.stdout
    out = run("check", "--sex", "female")
    assert out.returncode == 0 and "female_rule_coverage" in out.stdout and "modifier_data" in out.stdout
