import pytest

from humansim import Context, Registry, load_couplings, load_modifiers, load_sources, resolve, solve
from humansim.editing import (Effect, GeneEdit, catalogue_checks, edit_modifiers, load_catalogue)
from humansim.pharmacology import Dose, load_drugs, simulate
from humansim.registry import RegistryError


@pytest.fixture(scope="module")
def env():
    reg = Registry.load()
    return reg, load_modifiers(), load_couplings(reg)


def hyp(param="skeletal_muscle_mass", op="ratio", value=1.2):
    return Effect(param, op, value, tier="hypothetical", rationale="test")


def test_catalogue_is_consistent_and_every_effect_is_labelled(env):
    reg = env[0]
    cat = load_catalogue()
    assert catalogue_checks(reg, cat, load_sources()) == []
    assert cat  # starter set exists
    for loc in cat.values():
        for _, effects in loc.states.values():
            assert all(e.tier in ("measured", "hypothetical") for e in effects)
            assert all(e.rationale for e in effects if e.tier == "hypothetical")


def test_effect_validation():
    with pytest.raises(RegistryError):
        Effect("hemoglobin", "ratio", 1.1, tier="measured")           # measured needs a source
    with pytest.raises(RegistryError):
        Effect("hemoglobin", "ratio", 1.1, tier="hypothetical")       # hypothetical needs a rationale
    with pytest.raises(RegistryError):
        Effect("hemoglobin", "ratio", -1, tier="hypothetical", rationale="x")
    with pytest.raises(RegistryError):
        Effect("hemoglobin", "log", 1, tier="hypothetical", rationale="x")
    with pytest.raises(RegistryError):
        Effect("hemoglobin", "ratio", 1.1, tier="guess", rationale="x")


def test_edit_validation():
    with pytest.raises(ValueError):
        GeneEdit("X")                                                  # nothing to apply
    with pytest.raises(ValueError):
        GeneEdit("X", hypothesis=(hyp(),), mode="germline", fraction=0.5)
    with pytest.raises(ValueError):
        GeneEdit("X", hypothesis=(hyp(),), fraction=0)
    with pytest.raises(ValueError):
        GeneEdit("X", hypothesis=(hyp(),), mode="clone")
    with pytest.raises(ValueError):
        GeneEdit("X", locus="MSTN")


def test_hypothetical_edit_changes_value_and_is_flagged_everywhere(env):
    reg, mods, cps = env
    base = resolve(reg, Context(), mods)
    s = resolve(reg, Context(gene_edits=(GeneEdit("MSTN", hypothesis=(hyp(value=1.25),)),)), mods)
    assert s.values["skeletal_muscle_mass"] == pytest.approx(base.values["skeletal_muscle_mass"] * 1.25)
    f = s.flags["skeletal_muscle_mass"]
    assert {"hypothetical", "gene_edit", "modifier_hypothetical"} <= f
    assert "hypothetical" in s.flags["muscle_fraction"]              # inherited by the derived value
    assert "HYPOTHETICAL" in s.trace[-1].note


def test_somatic_fraction_scales_ratio_and_offset_and_is_flagged(env):
    reg, mods, cps = env
    base = resolve(reg, Context(), mods).values["skeletal_muscle_mass"]
    half = resolve(reg, Context(gene_edits=(GeneEdit("G", hypothesis=(hyp(value=1.4),), fraction=0.5),)), mods)
    assert half.values["skeletal_muscle_mass"] == pytest.approx(base * 1.2)
    assert "mosaic_assumed_linear" in half.flags["skeletal_muscle_mass"]
    off = edit_modifiers([GeneEdit("G", hypothesis=(hyp(op="offset", value=4.0),), fraction=0.25)])
    assert off[0].effects["skeletal_muscle_mass"] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        edit_modifiers([GeneEdit("G", hypothesis=(hyp(op="set", value=40.0),), fraction=0.5)])
    full = edit_modifiers([GeneEdit("G", hypothesis=(hyp(op="set", value=40.0),))])
    assert full[0].effects["skeletal_muscle_mass"] == 40.0


def test_germline_is_flagged_and_children_are_marked_unmodelled(env):
    reg, mods, cps = env
    e = GeneEdit("G", hypothesis=(hyp(),), mode="germline")
    adult = resolve(reg, Context(gene_edits=(e,)), mods)
    assert "germline" in adult.flags["skeletal_muscle_mass"]
    assert "age_effect_unmodelled" not in adult.flags["skeletal_muscle_mass"]
    child = resolve(reg, Context(age=10, gene_edits=(e,)), mods)
    assert "age_effect_unmodelled" in child.flags["skeletal_muscle_mass"]


def test_unknown_locus_state_and_param_are_refused(env):
    reg, mods, cps = env
    with pytest.raises(KeyError):
        resolve(reg, Context(gene_edits=(GeneEdit("X", locus="NOPE", state="a"),)), mods)
    with pytest.raises(KeyError):
        resolve(reg, Context(gene_edits=(GeneEdit("MSTN", locus="MSTN", state="nope"),)), mods)
    with pytest.raises(KeyError):
        resolve(reg, Context(gene_edits=(GeneEdit("X", hypothesis=(hyp(param="no_such_param"),)),)), mods)


def test_solver_propagates_a_hypothetical_edit_and_marks_downstream(env):
    reg, mods, cps = env
    base = solve(reg, Context(), mods, cps)
    s = solve(reg, Context(gene_edits=(GeneEdit("MSTN", locus="MSTN", state="loss_of_function"),)), mods, cps)
    assert s.solver["converged"]
    assert s.values["body_mass"] > base.values["body_mass"]
    assert "downstream_of_hypothetical" in s.flags["body_mass"]
    assert "downstream_of_hypothetical" in s.flags["body_mass_index"]    # derived from a solved value
    assert "downstream_of_hypothetical" not in base.flags["body_mass"]


def test_epor_edit_ripples_through_oxygen_content_and_cardiac_output(env):
    reg, mods, cps = env
    base = solve(reg, Context(), mods, cps)
    s = solve(reg, Context(gene_edits=(GeneEdit("EPOR", locus="EPOR", state="truncated_gain_of_function"),)), mods, cps)
    assert s.values["arterial_o2_content"] > base.values["arterial_o2_content"]
    assert s.values["cardiac_output"] < base.values["cardiac_output"]    # Fick: more O2 per litre, less flow needed


def test_edits_affect_drug_curves_and_are_flagged(env):
    reg, mods, cps = env
    drugs = load_drugs()
    e = GeneEdit("G", hypothesis=(hyp(param="liver_mass", value=1.5),))
    base = simulate(reg, solve(reg, Context(), mods, cps), drugs["acetaminophen"], [Dose(1000, 0, 0.25)], 24, mods)
    r = simulate(reg, solve(reg, Context(gene_edits=(e,)), mods, cps), drugs["acetaminophen"], [Dose(1000, 0, 0.25)], 24, mods)
    assert r.clearance_l_h > base.clearance_l_h
    assert "hypothetical_inputs" in r.flags and "hypothetical_inputs" not in base.flags


def test_reference_state_excludes_edits(env):
    reg, mods, cps = env
    e = GeneEdit("G", hypothesis=(hyp(),))
    s = solve(reg, Context(gene_edits=(e,)), mods, cps)
    assert "skeletal_muscle_mass" in {t.param for t in s.trace if t.stage == "genotype"}
    assert s.values["skeletal_muscle_mass"] == pytest.approx(29.0 * 1.2)  # couplings do not undo it


def test_gene_index_rejects_unknown_symbols_when_present():
    e = GeneEdit("NOTAGENE", hypothesis=(hyp(),))
    with pytest.raises(KeyError):
        edit_modifiers([e], gene_index={"MSTN": {}})
    assert edit_modifiers([e], gene_index={})  # no index imported: any symbol accepted
