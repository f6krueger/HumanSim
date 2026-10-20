import subprocess
import sys

from humansim import Context, resolve
from humansim.checks import range_flags, segment_stature_sum, state_checks, tissue_mass_budget


def test_tissue_budget_passes_for_baseline(registry):
    r = tissue_mass_budget(resolve(registry, Context()))
    assert r.status == "pass"
    assert 0 <= r.details["unassigned_kg"] < 3


def test_tissue_budget_fails_when_tissues_exceed_mass(registry):
    s = resolve(registry, Context(edits={"skeletal_muscle_mass": 100.0}))
    assert tissue_mass_budget(s).status == "fail"


def test_tissue_budget_warns_when_a_tissue_is_missing(registry):
    s = resolve(registry, Context())
    s.values["liver_mass"] = None  # simulate a tissue with no value
    assert tissue_mass_budget(s).status == "warn"


def test_segment_check_skipped_without_data_then_runs_with_data(registry):
    assert segment_stature_sum(resolve(registry, Context())).status == "skipped"
    good = Context(edits={"head_height": 23.0, "neck_length": 12.0, "torso_length": 55.0,
                          "crotch_height": 86.0})
    assert segment_stature_sum(resolve(registry, good)).status == "pass"
    bad = Context(edits={"head_height": 23.0, "neck_length": 12.0, "torso_length": 55.0,
                         "crotch_height": 60.0})
    assert segment_stature_sum(resolve(registry, bad)).status == "fail"


def test_range_flags_check(registry):
    assert range_flags(resolve(registry, Context())).status == "pass"
    assert range_flags(resolve(registry, Context(edits={"hemoglobin": 25.0}))).status == "warn"


def test_state_checks_run(registry):
    assert len(state_checks(resolve(registry, Context()))) == 3


def _cli(*args):
    return subprocess.run([sys.executable, "-m", "humansim", *args],
                          capture_output=True, text=True, timeout=60)


def test_cli_show_coverage_check():
    out = _cli("show")
    assert out.returncode == 0 and "body_mass_index" in out.stdout
    out = _cli("coverage", "--missing")
    assert out.returncode == 0 and "arm_span" in out.stdout
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
