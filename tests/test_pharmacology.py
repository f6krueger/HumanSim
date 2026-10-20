import math

import pytest

from humansim import Context, Registry, load_couplings, load_modifiers, load_sources, resolve, solve
from humansim.checks import drug_checks
from humansim.pharmacology import Dose, Drug, load_drugs, pk_parameters, simulate
from humansim.registry import RegistryError


@pytest.fixture(scope="module")
def reg():
    return Registry.load()


@pytest.fixture(scope="module")
def env(reg):
    return load_modifiers(), load_couplings(reg), load_drugs()


def test_drug_data_loads_with_sources(env):
    mods, cps, drugs = env
    assert set(drugs) == {"acetaminophen", "midazolam"}
    assert set(d.source_ref for d in drugs.values()) <= set(load_sources())
    assert all(r.status != "fail" for r in drug_checks(Registry.load(), drugs, load_sources(), mods))


def test_reference_man_uses_label_values_per_kg(reg, env):
    mods, cps, drugs = env
    s = resolve(reg, Context(), mods)
    cl, v, mass = pk_parameters(reg, s, drugs["acetaminophen"], mods)
    assert mass == pytest.approx(73.0)
    assert cl == pytest.approx(0.27 * 73.0) and v == pytest.approx(0.8 * 73.0)


def test_iv_bolus_matches_analytic_solution(reg, env):
    mods, cps, drugs = env
    s = resolve(reg, Context(), mods)
    r = simulate(reg, s, drugs["acetaminophen"], [Dose(1000.0, 0.0, 0.0)], 12.0, mods, dt=0.002)
    from humansim.dynamics import rk4  # independent numerical check of the closed form
    k0 = r.clearance_l_h / r.volume_l
    num = rk4(lambda t, y, u: {"A": (1000.0 / (1 / 60) if t < 1 / 60 else 0.0) - k0 * y["A"]}, {"A": 0.0}, 12.0, 1e-4)
    assert r.conc_mg_l[-1] == pytest.approx(num.final()["A"] / r.volume_l, rel=2e-2)
    k = r.clearance_l_h / r.volume_l
    t_end = 1 / 60
    c_end = 1000 / r.volume_l * (1 - math.exp(-k * t_end)) / (k * t_end)  # amount at end of the 1-min window
    i = min(range(len(r.times_h)), key=lambda j: abs(r.times_h[j] - 6.0))
    expect = c_end * math.exp(-k * (r.times_h[i] - t_end))
    assert r.conc_mg_l[i] == pytest.approx(expect, rel=1e-9)
    long = simulate(reg, s, drugs["acetaminophen"], [Dose(1000.0, 0.0, 0.0)], 48.0, mods)
    assert long.auc_mg_h_l() == pytest.approx(1000 / long.clearance_l_h, rel=0.005)  # AUC = dose / CL


def test_infusion_conserves_dose(reg, env):
    mods, cps, drugs = env
    s = resolve(reg, Context(), mods)
    r = simulate(reg, s, drugs["midazolam"], [Dose(10.0, 0.0, 2.0)], 60.0, mods)
    assert r.auc_mg_h_l() == pytest.approx(10.0 / r.clearance_l_h, rel=0.01)


def test_repeat_doses_accumulate(reg, env):
    mods, cps, drugs = env
    s = resolve(reg, Context(), mods)
    one = simulate(reg, s, drugs["acetaminophen"], [Dose(1000, 0, 0.25)], 30, mods)
    three = simulate(reg, s, drugs["acetaminophen"], [Dose(1000, t, 0.25) for t in (0, 6, 12)], 30, mods)
    assert three.cmax_mg_l > one.cmax_mg_l and three.auc_mg_h_l() == pytest.approx(3 * one.auc_mg_h_l(), rel=0.02)


def test_bigger_liver_clears_faster_and_is_flagged(reg, env):
    mods, cps, drugs = env
    base = solve(reg, Context(), mods, cps)
    big = solve(reg, Context(edits={"liver_mass": 2.5}), mods, cps)
    for d in drugs.values():
        a = simulate(reg, base, d, [Dose(100, 0, 0.1)], 24, mods)
        b = simulate(reg, big, d, [Dose(100, 0, 0.1)], 24, mods)
        assert b.clearance_l_h > a.clearance_l_h and b.half_life_h < a.half_life_h
        assert "scaled_from_changed_body" in b.flags and "scaled_from_changed_body" not in a.flags


def test_heavier_body_has_larger_volume_and_lower_peak(reg, env):
    mods, cps, drugs = env
    base = solve(reg, Context(), mods, cps)
    heavy = solve(reg, Context(edits={"body_mass": 120.0}), mods, cps)
    a = simulate(reg, base, drugs["acetaminophen"], [Dose(1000, 0, 0.25)], 12, mods)
    b = simulate(reg, heavy, drugs["acetaminophen"], [Dose(1000, 0, 0.25)], 12, mods)
    assert b.volume_l > a.volume_l and b.cmax_mg_l < a.cmax_mg_l


def test_women_get_smaller_absolute_clearance(reg, env):
    mods, cps, drugs = env
    m = pk_parameters(reg, resolve(reg, Context(), mods), drugs["midazolam"], mods)
    f = pk_parameters(reg, resolve(reg, Context(sex="female"), mods), drugs["midazolam"], mods)
    assert f[0] < m[0] and f[1] < m[1]


def test_unknown_route_and_bad_input(reg, env):
    mods, cps, drugs = env
    s = resolve(reg, Context(), mods)
    with pytest.raises(ValueError):
        simulate(reg, s, drugs["midazolam"], [], 10, mods)
    oral = Drug("x", "x", ("oral",), 1, 1, None, "unsourced", None, {}, {})
    with pytest.raises(ValueError):
        simulate(reg, s, oral, [Dose(1)], 10, mods)
    with pytest.raises(RegistryError):
        Drug("x", "x", ("iv",), 1, 1, None, "referenced", None, {}, {})
    with pytest.raises(RegistryError):
        Drug("x", "x", ("iv",), 1, 1, None, "unsourced", None, {"nonsense": 1.0}, {})
    with pytest.raises(RegistryError):
        Drug("x", "x", ("iv",), -1, 1, None, "unsourced", None, {}, {})
