import pytest
import yaml

from humansim import Registry, RegistryError
from humansim.validation import load_validation, validation_problems


@pytest.fixture(scope="module")
def female():
    return load_validation("ansur2_female")


def test_female_file_is_consistent_with_registry(registry, female):
    assert validation_problems(female, registry) == []
    assert female["sex"] == "female"
    assert female["source_ref"] == "ansur2"


def test_female_targets_are_not_baseline_values(registry):
    # Female data must stay a validation target; the sex rules derive the female body.
    assert registry["height"].baseline("female", "adult") is None


def test_every_validation_parameter_has_a_male_value(registry, female):
    for pid in female["parameters"]:
        assert registry[pid].baseline("male", "adult") is not None, pid


def test_ratio_to_male_matches_the_two_means(registry, female):
    for pid, rec in female["parameters"].items():
        male = registry[pid].baseline("male", "adult")
        assert rec["ratio_to_male"] == pytest.approx(rec["mean"] / male, abs=2e-3), pid


def test_sex_differences_have_the_expected_direction(female):
    p = female["parameters"]
    for pid in ("height", "body_mass", "arm_span", "hand_length", "foot_length", "biacromial_breadth"):
        assert p[pid]["ratio_to_male"] < 1, pid
    assert p["body_mass"]["ratio_to_male"] < p["height"]["ratio_to_male"]  # mass differs more than height


def test_male_values_are_in_plausible_ranges(registry):
    h = registry["height"]
    assert 170 < h.baseline("male", "adult") < 182
    assert 4 < h.spread("male", "adult") < 9
    lo, hi = h.ref_range["male"]
    flo, fhi = h.ref_range["female"]
    assert fhi < hi and flo < lo  # women are shorter
    assert lo < h.baseline("male", "adult") < hi
    assert h.status == "referenced" and h.source_ref == "ansur2"


def test_spread_accessor(registry):
    assert registry["body_mass"].spread("male", "adult") > 5
    assert registry["liver_mass"].spread("male", "adult") is None
    assert registry["height"].spread("female", "adult") is None


def test_loader_rejects_bad_sd(tmp_path):
    (tmp_path / "values").mkdir()

    def write(p):
        (tmp_path / "values" / "x.yaml").write_text(yaml.safe_dump({"group": "g", "parameters": [p]}))

    base = {"id": "a", "name": "A", "unit": "kg", "kind": "anchor", "values": {"male": {"adult": 1.0}}}
    write({**base, "sd": {"male": {"adult": -1.0}}})
    with pytest.raises(RegistryError, match="positive"):
        Registry.load(tmp_path)
    write({**base, "sd": {"female": {"adult": 1.0}}})
    with pytest.raises(RegistryError, match="no value"):
        Registry.load(tmp_path)


def test_problems_are_reported(registry, female):
    bad = {"parameters": {"not_a_param": {"unit": "cm", "mean": 1, "sd": 1, "p2_5": 0, "p97_5": 2, "n": 5},
                          "height": {"unit": "mm", "mean": 1620, "sd": 60, "p2_5": 1500, "p97_5": 1750, "n": 5},
                          "body_mass": {"unit": "kg", "mean": 99.0, "sd": 5, "p2_5": 50, "p97_5": 60, "n": 5}}}
    msgs = validation_problems(bad, registry)
    assert any("not a registry parameter" in m for m in msgs)
    assert any("unit" in m for m in msgs)
    assert any("outside its own" in m for m in msgs)
