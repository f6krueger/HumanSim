import pytest

from humansim import Context, Modifier, Parameter, Registry, load_modifiers, resolve
from humansim.checks import sex_rule_coverage, modifier_checks, range_flags
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


# --- both sexes are stored --------------------------------------------------------------------
def test_female_values_are_stored_not_derived(registry, female, female_doc):
    for pid, rec in female_doc["parameters"].items():
        assert registry[pid].baseline("female", "adult") == pytest.approx(rec["mean"], abs=0.006), pid
        assert female.get(pid) == registry[pid].baseline("female", "adult")
        assert "derived_from_male" not in female.flags[pid]
        assert registry[pid].spread("female", "adult") == pytest.approx(rec["sd"], abs=0.006)


def test_sex_rules_do_not_alter_stored_values(registry, mods):
    for sex in ("male", "female"):
        with_rules = resolve(registry, Context(sex=sex), mods)
        without = resolve(registry, Context(sex=sex), [])
        assert with_rules.get("height") == without.get("height")
        assert with_rules.get("body_mass") == without.get("body_mass")
    assert resolve(registry, Context(sex="male"), mods).trace == []


def test_rules_agree_with_stored_pairs(registry, mods):
    assert modifier_checks(registry, mods)[0].status == "pass"
    bad = Modifier("x", "sex", "ratio", {"height": 0.5}, when={"sex": "female"})
    assert modifier_checks(registry, [bad])[0].status == "fail"


def _only_one_sex(sex):
    a = Parameter(id="a", name="A", unit="kg", kind="anchor", group="g", status="referenced",
                  source_ref="s", values={sex: {"adult": 100.0 if sex == "male" else 60.0}})
    return Registry({"a": a})


def test_female_derived_from_male_and_male_derived_from_female_are_inverse():
    rule = Modifier("r", "sex", "ratio", {"a": 0.6}, when={"sex": "female"})
    f = resolve(_only_one_sex("male"), Context(sex="female"), [rule])
    assert f.get("a") == pytest.approx(60.0) and "derived_from_male" in f.flags["a"]
    m = resolve(_only_one_sex("female"), Context(sex="male"), [rule])
    assert m.get("a") == pytest.approx(100.0) and "derived_from_female" in m.flags["a"]
    assert "inverse" in m.trace[0].note
    # round trip
    assert m.get("a") * 0.6 == pytest.approx(60.0)


def test_female_only_value_without_rule_is_empty_for_male():
    s = resolve(_only_one_sex("female"), Context(sex="male"), [])
    assert s.get("a") is None and s.flags["a"] == {"missing", "no_sex_rule"}
    kept = resolve(_only_one_sex("female"), Context(sex="male", keep_unadjusted=True), [])
    assert kept.get("a") == 60.0 and "derived_from_female" in kept.flags["a"]


def test_set_rule_cannot_be_inverted_for_male():
    rule = Modifier("r", "sex", "set", {"a": 60.0}, when={"sex": "female"})
    s = resolve(_only_one_sex("female"), Context(sex="male"), [rule])
    assert s.get("a") is None and "no_sex_rule" in s.flags["a"]


def test_every_sex_rule_targets_a_valid_parameter_once(registry, mods):
    assert modifier_checks(registry, mods)[0].status == "pass"
    assert len(mods) > 150  # one rule per stored male/female pair, plus assumed-neutral ones
    assert {m.stage for m in mods} == {"sex", "hormone"}
    assert all(m.stage == "hormone" for m in mods if m.rule_type == "hormonal")
    assert all(m.when == {"sex": "female"} for m in mods if m.stage == "sex")


def test_rule_types_are_labelled_and_sourced(mods):
    types = {m.rule_type for m in mods}
    assert types == {"scaling", "anatomical", "physiological", "assumed_neutral", "hormonal"}
    for m in mods:
        assert m.status in ("referenced", "unsourced")
        if m.status == "referenced":
            assert m.source_ref
        if m.rule_type == "assumed_neutral":
            assert m.status == "unsourced"


# --- the ICRP female/male pairs -------------------------------------------------------------
def test_icrp_female_values(female, registry):
    assert registry["hematocrit"].baseline("female", "adult") == 0.42
    assert female.get("hematocrit") == pytest.approx(0.42, abs=1e-3)
    assert female.get("cardiac_output") == pytest.approx(5.9, abs=1e-3)
    assert female.get("blood_volume") == pytest.approx(3.9, abs=1e-3)


# --- parameters without a rule are not silently filled with male values ---------------------
def test_no_rule_means_no_value_for_female(female, male):
    assert male.get("brain_metabolic_load") is not None
    assert female.get("brain_metabolic_load") is None
    assert female.flags["brain_metabolic_load"] == {"missing", "no_sex_rule"}


def test_keep_unadjusted_keeps_the_male_value_but_flags_it(registry, mods, male):
    s = resolve(registry, Context(sex="female", keep_unadjusted=True), mods)
    assert s.get("brain_metabolic_load") == male.get("brain_metabolic_load")
    assert {"derived_from_male", "no_sex_rule"} <= s.flags["brain_metabolic_load"]


def test_female_hemoglobin_follows_hematocrit_ratio(female, male):
    assert female.get("hemoglobin") == pytest.approx(male.get("hemoglobin") * 0.42 / 0.45, rel=1e-3)
    assert 12.0 < female.get("hemoglobin") < 16.0
    assert "modifier_unsourced" in female.flags["hemoglobin"]


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
    s = resolve(registry, Context(sex="female", edits={"body_mass": 80.0, "brain_metabolic_load": 15.0}), mods)
    assert s.get("body_mass") == 80.0
    assert s.get("brain_metabolic_load") == 15.0 and s.flags["brain_metabolic_load"] == {"edited"}
    assert [t.stage for t in s.trace if t.param == "body_mass"] == ["sex", "edit"]


def test_female_coverage_check_lists_missing_rules(female):
    r = sex_rule_coverage(female)
    assert r.status == "warn" and "brain_metabolic_load" in r.details["no_sex_rule"]
    male_cov = sex_rule_coverage(resolve(Registry.load(), Context()))
    assert "prolactin" in male_cov.details["no_sex_rule"]  # a female-only hormone has no male value


# --- reference ranges are per sex -------------------------------------------------------------
def test_female_baseline_is_inside_female_ranges_not_flagged(female):
    assert range_flags(female).status == "pass"
    assert "outside_reference_range" not in female.flags["neck_circumference"]


def test_male_range_is_not_applied_to_a_female(registry, mods):
    # 185 cm is inside the male range (up to 189.7) but far above the female 97.5th percentile (175.8).
    assert "outside_reference_range" not in resolve(registry, Context(edits={"height": 185.0}), mods).flags["height"]
    assert "outside_reference_range" in resolve(registry, Context(sex="female", edits={"height": 185.0}), mods).flags["height"]
    # Haemoglobin: 14 g/dL is inside the male range (13-18) but 18 is outside the female range (11.5-16.5).
    assert "outside_reference_range" not in resolve(registry, Context(edits={"hemoglobin": 17.5}), mods).flags["hemoglobin"]
    assert "outside_reference_range" in resolve(registry, Context(sex="female", edits={"hemoglobin": 17.5}), mods).flags["hemoglobin"]


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
    assert out.returncode == 0 and "Independent check" in out.stdout and "body_mass" in out.stdout
    out = run("show", "--sex", "female", "--group", "anthropometry", "--trace")
    assert out.returncode == 0 and "sex_female_height" in out.stdout
    out = run("check", "--sex", "female")
    assert out.returncode == 0 and "sex_rule_coverage" in out.stdout and "modifier_data" in out.stdout
