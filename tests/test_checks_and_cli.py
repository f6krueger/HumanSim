import subprocess
import sys

from humansim import Context, resolve
from humansim.checks import range_flags, segment_stature_sum, state_checks, tissue_mass_budget


def test_tissue_budget_warns_about_population_mismatch_at_baseline(registry):
    r = tissue_mass_budget(resolve(registry, Context()))
    assert r.status == "warn"
    assert "different populations" in r.message
    assert r.details["missing"] == []
    assert r.details["unassigned_kg"] > 8


def test_tissue_budget_passes_when_mass_matches_the_organ_data(registry):
    r = tissue_mass_budget(resolve(registry, Context(edits={"body_mass": 73.0})))
    assert r.status == "pass"
    assert 0 <= r.details["unassigned_kg"] < 3


def test_tissue_budget_fails_when_tissues_exceed_mass(registry):
    s = resolve(registry, Context(edits={"skeletal_muscle_mass": 100.0}))
    assert tissue_mass_budget(s).status == "fail"


def test_tissue_budget_reports_a_missing_tissue(registry):
    s = resolve(registry, Context())
    s.values["liver_mass"] = None  # simulate a tissue with no value
    r = tissue_mass_budget(s)
    assert r.status == "warn" and r.details["missing"] == ["liver_mass"]


def test_segment_check_skipped_without_data_then_runs_with_data(registry):
    base = resolve(registry, Context())
    assert segment_stature_sum(base).status == "skipped"  # head_height, neck_length have no value
    need = base.get("height") - base.get("torso_length") - base.get("crotch_height")
    good = Context(edits={"head_height": need - 9.0, "neck_length": 9.0})
    assert segment_stature_sum(resolve(registry, good)).status == "pass"
    bad = Context(edits={"head_height": need - 9.0, "neck_length": 9.0, "crotch_height": 60.0})
    assert segment_stature_sum(resolve(registry, bad)).status == "fail"


def test_range_flags_check(registry):
    assert range_flags(resolve(registry, Context())).status == "pass"
    assert range_flags(resolve(registry, Context(edits={"hemoglobin": 25.0}))).status == "warn"
    assert range_flags(resolve(registry, Context(edits={"height": 215.0}))).status == "warn"


def test_state_checks_run(registry):
    assert len(state_checks(resolve(registry, Context()))) == 6


def _cli(*args):
    return subprocess.run([sys.executable, "-m", "humansim", *args],
                          capture_output=True, text=True, timeout=60)


def test_cli_show_coverage_check():
    out = _cli("show")
    assert out.returncode == 0 and "body_mass_index" in out.stdout
    out = _cli("coverage", "--missing")
    assert out.returncode == 0 and "femur_length" in out.stdout
    out = _cli("check")
    assert out.returncode == 0 and "tissue_mass_budget" in out.stdout


def test_cli_edit_with_units_and_changed_view():
    out = _cli("show", "--set", "body_mass=198 lb", "--changed")
    assert out.returncode == 0
    assert "body_mass" in out.stdout and "body_mass_index" in out.stdout
    assert "height" not in out.stdout.split()  # unchanged values are hidden


def test_cli_bad_input_exits_cleanly():
    out = _cli("show", "--set", "nonsense=1")
    assert out.returncode != 0 and "Traceback" not in out.stderr
    out = _cli("show", "--set", "body_mass_index=30")
    assert out.returncode != 0 and "derived" in out.stderr and "Traceback" not in out.stderr
