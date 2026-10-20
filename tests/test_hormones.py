import math

import pytest

from humansim import Context, Registry, load_couplings, load_modifiers, resolve, solve
from humansim import cycle, hpg
from humansim.dynamics import rk4


@pytest.fixture(scope="module")
def reg():
    return Registry.load()


@pytest.fixture(scope="module")
def mods():
    return load_modifiers()


@pytest.fixture(scope="module")
def cps(reg):
    return load_couplings(reg)


def test_hormone_data_and_phase_lookup(reg, mods):
    assert reg["testosterone_total"].baseline("male", "adult") == pytest.approx(477.49, abs=0.01)
    luteal = resolve(reg, Context(sex="female", cycle_phase="luteal"), mods)
    follicular = resolve(reg, Context(sex="female"), mods)
    assert luteal.values["progesterone"] > 10 * follicular.values["progesterone"]
    assert luteal.phase == "luteal"
    # a luteal progesterone is normal for the luteal phase but would be high for the follicular phase
    assert "outside_reference_range" not in luteal.flags["progesterone"]
    high = resolve(reg, Context(sex="female", cycle_phase="follicular", hormone_state={"progesterone": 9.24}), mods)
    assert "outside_reference_range" in high.flags["progesterone"]


def test_hormone_state_is_flagged_and_range_checked(reg, mods):
    s = resolve(reg, Context(hormone_state={"testosterone_total": 1345}), mods)
    assert s.values["testosterone_total"] == 1345
    assert {"hormone_set", "outside_reference_range"} <= s.flags["testosterone_total"]
    with pytest.raises(KeyError):
        resolve(reg, Context(hormone_state={"body_mass": 80}), mods)
    with pytest.raises(ValueError):
        resolve(reg, Context(hormone_state={"testosterone_total": -1}), mods)


def test_testosterone_effect_on_muscle_follows_log_slope(reg, mods):
    base = resolve(reg, Context(), mods)
    ref = base.values["testosterone_total"]
    s = resolve(reg, Context(hormone_state={"testosterone_total": ref * math.e}), mods)
    assert s.values["skeletal_muscle_mass"] - base.values["skeletal_muscle_mass"] == pytest.approx(2.9499, abs=1e-3)
    assert "hormone_effect" in s.flags["skeletal_muscle_mass"]


def test_no_hormone_change_means_no_effect(reg, mods):
    base = resolve(reg, Context(), mods)
    same = resolve(reg, Context(hormone_state={"testosterone_total": base.values["testosterone_total"]}), mods)
    assert same.values["skeletal_muscle_mass"] == pytest.approx(base.values["skeletal_muscle_mass"])


def test_effect_outside_studied_concentrations_is_flagged(reg, mods):
    inside = resolve(reg, Context(hormone_state={"testosterone_total": 1345}), mods)
    outside = resolve(reg, Context(hormone_state={"testosterone_total": 5000}), mods)
    assert "outside_validity" not in inside.flags["skeletal_muscle_mass"]
    assert any("outside" in f for f in outside.flags["skeletal_muscle_mass"])


def test_effect_does_not_apply_to_women(reg, mods):
    f = resolve(reg, Context(sex="female"), mods)
    f2 = resolve(reg, Context(sex="female", hormone_state={"testosterone_total": 60}), mods)
    assert f2.values["skeletal_muscle_mass"] == f.values["skeletal_muscle_mass"]


def test_solver_carries_hormone_effect_to_dependents(reg, mods, cps):
    base = solve(reg, Context(), mods, cps)
    s = solve(reg, Context(hormone_state={"testosterone_total": 1345}), mods, cps)
    assert s.solver["converged"]
    assert s.values["skeletal_muscle_mass"] > base.values["skeletal_muscle_mass"]  # not undone by couplings
    assert s.values["body_mass"] - base.values["body_mass"] == pytest.approx(
        s.values["skeletal_muscle_mass"] - base.values["skeletal_muscle_mass"], rel=0.05)
    assert s.values["resting_metabolic_rate"] > base.values["resting_metabolic_rate"]


def test_solver_unedited_state_is_still_the_reference(reg, mods, cps):
    a, b = solve(reg, Context(), mods, cps), resolve(reg, Context(), mods)
    for k, v in b.values.items():
        assert a.values[k] == pytest.approx(v)


def test_rk4_matches_exponential_decay():
    tr = rk4(lambda t, y, u: {"x": -y["x"]}, {"x": 1.0}, 2.0, 0.01)
    assert tr.final()["x"] == pytest.approx(math.exp(-2), rel=1e-8)
    assert tr.t[-1] == pytest.approx(2.0)


def test_rk4_flags_blow_up():
    tr = rk4(lambda t, y, u: {"x": y["x"] ** 3}, {"x": 10.0}, 5.0, 0.01)
    assert not tr.converged


def test_hpg_stays_at_reference_without_input(reg):
    r = hpg.simulate(reg, 100)
    assert r.final()["T"] == pytest.approx(r.reference["T"], rel=1e-9)
    assert r.final()["LH"] == pytest.approx(r.reference["LH"], rel=1e-9)
    assert r.assumed


def test_hpg_recovers_from_perturbation(reg):
    r = hpg.simulate(reg, 300, start={"T": 1345.0, "LH": 0.5})
    assert r.final()["T"] == pytest.approx(r.reference["T"], rel=1e-3)
    assert r.final()["LH"] == pytest.approx(r.reference["LH"], rel=1e-3)


def test_hpg_dosing_suppresses_lh_and_recovers(reg):
    r = hpg.simulate(reg, 400, dose_rate=lambda t: 20000.0 if t < 150 else 0.0)
    i = round(100 / 0.01)
    assert r.trajectory.y[i]["LH"] < 0.3 * r.reference["LH"]
    assert r.trajectory.y[i]["C"] < 0.5
    assert r.final()["T"] == pytest.approx(r.reference["T"], rel=1e-3)


def test_hpg_suppression_grows_with_dose(reg):
    lhs = [hpg.simulate(reg, 60, dose_rate=d).final()["LH"] for d in (0.0, 5000.0, 10000.0, 20000.0)]
    assert lhs == sorted(lhs, reverse=True) and lhs[0] > lhs[-1]


def test_cycle_profile_hits_phase_values_and_wraps(reg):
    for ph, day in cycle.ANCHOR_DAYS.items():
        assert cycle.level(reg, "lh", day) == pytest.approx(reg["lh"].baseline("female", f"adult:{ph}"))
    assert cycle.level(reg, "lh", 1) == pytest.approx(cycle.level(reg, "lh", 29))
    prof = cycle.profile(reg)
    assert max(prof, key=lambda d: prof[d]["lh"]) == 14.0
    assert max(prof, key=lambda d: prof[d]["progesterone"]) == 21.0


def test_cycle_context_sets_hormones_and_phase(reg, mods):
    ctx = cycle.cycle_context(reg, 21)
    s = resolve(reg, ctx, mods)
    assert s.phase == "luteal" and ctx.sex == "female"
    assert s.values["progesterone"] == pytest.approx(9.2396)
    assert "interpolated" in s.flags["progesterone"]
    assert cycle.phase_of(14) == "midcycle" and cycle.phase_of(3) == "follicular"
