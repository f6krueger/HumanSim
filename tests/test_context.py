import pytest

from humansim import Context, Modifier, resolve


def _base(registry):
    return resolve(registry, Context())


def test_baseline_male(registry):
    s = _base(registry)
    h, m = s.get("height"), s.get("body_mass")
    assert 170 < h < 182 and 75 < m < 95
    assert s.get("body_mass_index") == pytest.approx(m / (h / 100) ** 2)
    assert "missing" in s.flags["femur_length"]
    assert "unsourced" in s.flags["liver_mass"]
    assert "missing" in s.flags["intermembral_index"]
    assert s.missing_inputs["intermembral_index"] == ["femur_length", "tibia_length"]
    assert s.get("ape_index") == pytest.approx(s.get("arm_span") / h)


def test_edit_propagates_to_derived(registry):
    base = _base(registry)
    h = base.get("height")
    s = resolve(registry, Context(edits={"body_mass": 90.0}))
    assert s.get("body_mass_index") == pytest.approx(90 / (h / 100) ** 2)
    assert s.get("body_surface_area") > base.get("body_surface_area")
    changed = s.changed_from(base)
    assert "body_mass_index" in changed and "body_mass" in changed
    assert "height" not in changed
    assert "edited" in s.flags["body_mass"]
    assert "edited" in s.flags["body_mass_index"]  # derived inherits


def test_edit_fills_missing_values(registry):
    s = resolve(registry, Context(edits={"femur_length": 45.0, "tibia_length": 38.0}))
    expected = 100 * (s.get("upper_arm_length") + s.get("forearm_length")) / (45.0 + 38.0)
    assert s.get("intermembral_index") == pytest.approx(expected)
    assert s.flags["femur_length"] == {"edited"}
    assert "edited" in s.flags["intermembral_index"]


def test_cannot_edit_derived_or_unknown(registry):
    with pytest.raises(ValueError, match="derived"):
        resolve(registry, Context(edits={"body_mass_index": 30.0}))
    with pytest.raises(KeyError):
        resolve(registry, Context(edits={"no_such_parameter": 1.0}))
    with pytest.raises(ValueError, match="finite"):
        resolve(registry, Context(edits={"body_mass": float("nan")}))


def test_unsupported_sex(registry):
    with pytest.raises(ValueError):
        resolve(registry, Context(sex="other"))


def test_reference_range_flags(registry):
    s = resolve(registry, Context(edits={"hematocrit": 0.60}))
    assert "outside_reference_range" in s.flags["hematocrit"]
    tall = resolve(registry, Context(edits={"height": 215.0}))
    assert "outside_reference_range" in tall.flags["height"]  # far above the ANSUR II 97.5th percentile
    ok = _base(registry)
    assert "outside_reference_range" not in ok.flags["hematocrit"]
    assert "outside_reference_range" not in ok.flags["height"]


def test_edits_do_not_leak_between_calls(registry):
    before = _base(registry).get("body_mass")
    resolve(registry, Context(edits={"body_mass": 200.0}))
    assert _base(registry).get("body_mass") == before


# --- modifiers -------------------------------------------------------------------------
def test_modifier_ops(registry):
    base = _base(registry)
    mods = [
        Modifier("r", "sex", "ratio", {"body_mass": 2.0}),
        Modifier("o", "drug", "offset", {"height": 4.0}),
        Modifier("s", "genotype", "set", {"core_temperature": 38.0}),
    ]
    s = resolve(registry, Context(), mods)
    assert s.get("body_mass") == pytest.approx(2 * base.get("body_mass"))
    assert s.get("height") == pytest.approx(base.get("height") + 4.0)
    assert s.get("core_temperature") == 38.0
    assert {"modified"} <= s.flags["body_mass"]
    assert "modifier_unsourced" in s.flags["body_mass"]


def test_stage_order_is_fixed_not_list_order(registry):
    # sex ratio x2 must run before drug offset +1, whatever order they are listed in.
    base = _base(registry).get("body_mass")
    a = Modifier("a", "drug", "offset", {"body_mass": 1.0})
    b = Modifier("b", "sex", "ratio", {"body_mass": 2.0})
    assert resolve(registry, Context(), [a, b]).get("body_mass") == pytest.approx(2 * base + 1)
    assert resolve(registry, Context(), [b, a]).get("body_mass") == pytest.approx(2 * base + 1)


def test_edits_apply_after_modifiers(registry):
    m = Modifier("m", "sex", "ratio", {"body_mass": 2.0})
    s = resolve(registry, Context(edits={"body_mass": 80.0}), [m])
    assert s.get("body_mass") == 80.0
    assert [t.stage for t in s.trace] == ["sex", "edit"]


def test_when_conditions(registry):
    base = _base(registry).get("body_mass")
    f = Modifier("f", "sex", "ratio", {"body_mass": 0.5}, when={"sex": "female"})
    assert resolve(registry, Context(sex="male"), [f]).get("body_mass") == base
    g = Modifier("g", "genotype", "ratio", {"body_mass": 0.5}, when={"genotype:CYP2D6": "poor"})
    assert resolve(registry, Context(genotype={"CYP2D6": "poor"}), [g]).get("body_mass") == pytest.approx(base / 2)
    assert resolve(registry, Context(genotype={"CYP2D6": "normal"}), [g]).get("body_mass") == base
    d = Modifier("d", "drug", "offset", {"core_temperature": -1.0}, when={"drug:cooling": True})
    assert resolve(registry, Context(drugs={"cooling": 1.0}), [d]).get("core_temperature") == 36.0
    assert resolve(registry, Context(), [d]).get("core_temperature") == 37.0


def test_modifier_skips_missing_baseline_and_says_so(registry):
    m = Modifier("m", "sex", "ratio", {"femur_length": 1.1})
    s = resolve(registry, Context(), [m])
    assert s.get("femur_length") is None
    assert "modifier_skipped_no_baseline" in s.flags["femur_length"]


def test_modifier_on_derived_or_unknown_is_an_error(registry):
    with pytest.raises(KeyError):
        resolve(registry, Context(), [Modifier("m", "sex", "ratio", {"body_mass_index": 1.1})])


def test_modifier_validation():
    from humansim import RegistryError
    with pytest.raises(RegistryError):
        Modifier("m", "nonsense", "ratio", {"body_mass": 1.0})
    with pytest.raises(RegistryError):
        Modifier("m", "sex", "multiply", {"body_mass": 1.0})
    with pytest.raises(RegistryError):
        Modifier("m", "sex", "ratio", {})
