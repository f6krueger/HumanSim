import pytest

from humansim import Context, Modifier, resolve


def test_baseline_male(registry):
    s = resolve(registry, Context())
    assert s.get("height") == 176.0
    assert s.get("body_mass_index") == pytest.approx(73 / 1.76**2)
    assert "missing" in s.flags["arm_span"]
    assert "unsourced" in s.flags["liver_mass"]
    assert "missing" in s.flags["ape_index"]
    assert s.missing_inputs["ape_index"] == ["arm_span"]


def test_edit_propagates_to_derived(registry):
    base = resolve(registry, Context())
    s = resolve(registry, Context(edits={"body_mass": 90.0}))
    assert s.get("body_mass_index") == pytest.approx(90 / 1.76**2)
    assert s.get("body_surface_area") > base.get("body_surface_area")
    changed = s.changed_from(base)
    assert "body_mass_index" in changed and "body_mass" in changed
    assert "height" not in changed
    assert "edited" in s.flags["body_mass"]
    assert "edited" in s.flags["body_mass_index"]  # derived inherits


def test_edit_fills_missing_value(registry):
    s = resolve(registry, Context(edits={"arm_span": 180.0}))
    assert s.get("ape_index") == pytest.approx(180 / 176)
    assert s.flags["arm_span"] == {"edited"}


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


def test_reference_range_flag(registry):
    s = resolve(registry, Context(edits={"hematocrit": 0.60}))
    assert "outside_reference_range" in s.flags["hematocrit"]
    ok = resolve(registry, Context())
    assert "outside_reference_range" not in ok.flags["hematocrit"]


def test_edits_do_not_leak_between_calls(registry):
    resolve(registry, Context(edits={"body_mass": 200.0}))
    assert resolve(registry, Context()).get("body_mass") == 73.0


# --- modifiers -------------------------------------------------------------------------
def test_modifier_ops(registry):
    mods = [
        Modifier("r", "sex", "ratio", {"body_mass": 2.0}),
        Modifier("o", "drug", "offset", {"height": 4.0}),
        Modifier("s", "genotype", "set", {"core_temperature": 38.0}),
    ]
    s = resolve(registry, Context(), mods)
    assert s.get("body_mass") == 146.0
    assert s.get("height") == 180.0
    assert s.get("core_temperature") == 38.0
    assert {"modified"} <= s.flags["body_mass"]
    assert "modifier_unsourced" in s.flags["body_mass"]


def test_stage_order_is_fixed_not_list_order(registry):
    # sex ratio x2 must run before drug offset +1, whatever order they are listed in.
    a = Modifier("a", "drug", "offset", {"body_mass": 1.0})
    b = Modifier("b", "sex", "ratio", {"body_mass": 2.0})
    assert resolve(registry, Context(), [a, b]).get("body_mass") == 147.0
    assert resolve(registry, Context(), [b, a]).get("body_mass") == 147.0


def test_edits_apply_after_modifiers(registry):
    m = Modifier("m", "sex", "ratio", {"body_mass": 2.0})
    s = resolve(registry, Context(edits={"body_mass": 80.0}), [m])
    assert s.get("body_mass") == 80.0
    assert [t.stage for t in s.trace] == ["sex", "edit"]


def test_when_conditions(registry):
    f = Modifier("f", "sex", "ratio", {"body_mass": 0.5}, when={"sex": "female"})
    assert resolve(registry, Context(sex="male"), [f]).get("body_mass") == 73.0
    g = Modifier("g", "genotype", "ratio", {"body_mass": 0.5}, when={"genotype:CYP2D6": "poor"})
    assert resolve(registry, Context(genotype={"CYP2D6": "poor"}), [g]).get("body_mass") == 36.5
    assert resolve(registry, Context(genotype={"CYP2D6": "normal"}), [g]).get("body_mass") == 73.0
    d = Modifier("d", "drug", "offset", {"core_temperature": -1.0}, when={"drug:cooling": True})
    assert resolve(registry, Context(drugs={"cooling": 1.0}), [d]).get("core_temperature") == 36.0
    assert resolve(registry, Context(), [d]).get("core_temperature") == 37.0


def test_modifier_skips_missing_baseline_and_says_so(registry):
    m = Modifier("m", "sex", "ratio", {"arm_span": 1.1})
    s = resolve(registry, Context(), [m])
    assert s.get("arm_span") is None
    assert "modifier_skipped_no_baseline" in s.flags["arm_span"]


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
