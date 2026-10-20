"""Build the anthropometry data from the public ANSUR II CSV files.

Usage (run from anywhere; the CSVs are treated as untrusted input, so use -I):
    python -I tools/build_ansur_data.py ANSUR_II_MALE_Public.csv ANSUR_II_FEMALE_Public.csv

Writes:
    humansim/data/values/anthropometry.yaml      male mean, SD and 2.5-97.5 percentile range
    humansim/data/validation/ansur2_female.yaml  female statistics, the target that the Phase 2
                                                 sex rules must reproduce (not used as baseline)

Only the standard library is used. ANSUR II stores lengths in millimetres and body weight in
tenths of a kilogram; both are converted here. Nothing is typed in by hand: every number in the
output is computed from the files.
"""

from __future__ import annotations

import csv
import statistics
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_VALUES = ROOT / "humansim" / "data" / "values" / "anthropometry.yaml"
OUT_FEMALE = ROOT / "humansim" / "data" / "validation" / "ansur2_female.yaml"

MM_TO_CM = 0.1
TENTH_KG = 0.1

# (id, display name, unit, spec, note)
# spec: None                       -> no ANSUR II equivalent, parameter stays without a value
#       ("col", column, scale)     -> one column times scale
#       ("diff", col_a, col_b)     -> per-subject (col_a - col_b) in cm, then summarised
P = [
    ("height", "Height", "cm", ("col", "stature", MM_TO_CM), ""),
    ("body_mass", "Body mass", "kg", ("col", "weightkg", TENTH_KG), ""),
    ("body_density", "Whole-body density", "g/cm^3", None, "ANSUR II has no densitometry."),
    ("leg_length_trochanteric", "Leg length, trochanteric", "cm", ("col", "trochanterionheight", MM_TO_CM), "Height of the trochanterion above the floor."),
    ("crotch_height", "Crotch height", "cm", ("col", "crotchheight", MM_TO_CM), ""),
    ("hip_joint_height", "Hip joint centre height", "cm", None, "No direct ANSUR II measurement (trochanterion height is close but not the same point)."),
    ("sitting_height", "Sitting height, crown to rump", "cm", ("col", "sittingheight", MM_TO_CM), ""),
    ("crown_to_crotch", "Crown to crotch", "cm", ("diff", "stature", "crotchheight"), "Computed per subject as stature minus crotch height."),
    ("arm_span", "Arm span", "cm", ("col", "span", MM_TO_CM), ""),
    ("arm_length", "Arm length, acromion to fingertip", "cm", None, "No direct ANSUR II measurement."),
    ("standing_reach", "Standing reach", "cm", None, "ANSUR II only measures overhead reach while sitting."),
    ("hand_length", "Hand length", "cm", ("col", "handlength", MM_TO_CM), ""),
    ("foot_length", "Foot length", "cm", ("col", "footlength", MM_TO_CM), ""),
    ("biacromial_breadth", "Biacromial breadth", "cm", ("col", "biacromialbreadth", MM_TO_CM), ""),
    ("bi_iliac_breadth", "Bi-iliac breadth", "cm", ("col", "bicristalbreadth", MM_TO_CM), "ANSUR II bicristal breadth."),
    ("chest_circumference", "Chest circumference", "cm", ("col", "chestcircumference", MM_TO_CM), ""),
    ("chest_depth", "Chest depth", "cm", ("col", "chestdepth", MM_TO_CM), ""),
    ("waist_circumference", "Waist circumference", "cm", ("col", "waistcircumference", MM_TO_CM), ""),
    ("hip_circumference", "Hip circumference", "cm", ("col", "buttockcircumference", MM_TO_CM), "ANSUR II buttock circumference, used as hip circumference."),
    ("neck_circumference", "Neck circumference", "cm", ("col", "neckcircumference", MM_TO_CM), ""),
    ("upper_arm_circ_relaxed", "Upper arm circumference, relaxed", "cm", None, "ANSUR II only measures the flexed biceps."),
    ("upper_arm_circ_flexed", "Upper arm circumference, flexed", "cm", ("col", "bicepscircumferenceflexed", MM_TO_CM), ""),
    ("forearm_circumference", "Forearm circumference", "cm", ("col", "forearmcircumferenceflexed", MM_TO_CM), "ANSUR II flexed forearm circumference (maximum), not strictly the proximal site."),
    ("wrist_circumference", "Wrist circumference", "cm", ("col", "wristcircumference", MM_TO_CM), ""),
    ("thigh_circumference", "Thigh circumference, proximal", "cm", ("col", "thighcircumference", MM_TO_CM), ""),
    ("knee_circumference", "Knee circumference", "cm", None, "No ANSUR II knee circumference."),
    ("calf_circumference", "Calf circumference", "cm", ("col", "calfcircumference", MM_TO_CM), ""),
    ("ankle_circumference", "Ankle circumference", "cm", ("col", "anklecircumference", MM_TO_CM), ""),
    ("penis_length_erect", "Penis length, erect", "cm", None, "Not measured by ANSUR II."),
    ("penis_girth_erect", "Penis girth, erect", "cm", None, "Not measured by ANSUR II."),
    ("head_circumference", "Head circumference", "cm", ("col", "headcircumference", MM_TO_CM), ""),
    ("head_height", "Head height", "cm", None, "ANSUR II measures tragion-to-vertex, not full head height."),
    ("head_length", "Head length", "cm", ("col", "headlength", MM_TO_CM), ""),
    ("head_breadth", "Head breadth", "cm", ("col", "headbreadth", MM_TO_CM), ""),
    ("neck_length", "Neck length", "cm", None, "No ANSUR II neck length."),
    ("shoulder_height", "Shoulder height", "cm", ("col", "acromialheight", MM_TO_CM), "ANSUR II acromial height."),
    ("torso_length", "Torso length, suprasternal to crotch", "cm", ("diff", "suprasternaleheight", "crotchheight"), "Computed per subject as suprasternale height minus crotch height."),
    ("femur_length", "Femur length", "cm", None, "No ANSUR II femur length (bone lengths need imaging)."),
    ("tibia_length", "Tibia length", "cm", None, "No ANSUR II tibia length (bone lengths need imaging)."),
    ("foot_breadth", "Foot breadth", "cm", ("col", "footbreadthhorizontal", MM_TO_CM), ""),
    ("ankle_height", "Ankle height", "cm", ("col", "lateralmalleolusheight", MM_TO_CM), "Height of the lateral malleolus."),
    ("upper_arm_length", "Upper arm length", "cm", ("col", "acromionradialelength", MM_TO_CM), "Acromion to radiale (a surface landmark length, not humerus length)."),
    ("forearm_length", "Forearm length", "cm", ("col", "radialestylionlength", MM_TO_CM), "Radiale to stylion (a surface landmark length, not radius length)."),
    ("hand_breadth", "Hand breadth", "cm", ("col", "handbreadth", MM_TO_CM), ""),
]


def read(path: str) -> list[dict]:
    with open(path, newline="", encoding="latin-1") as fh:
        return list(csv.DictReader(fh))


def column(rows, name: str, scale: float) -> list[float]:
    if name not in rows[0]:
        raise SystemExit(f"Column {name!r} not found in the CSV")
    out = []
    for r in rows:
        txt = (r[name] or "").strip()
        if txt:
            out.append(float(txt) * scale)
    return out


def diff(rows, a: str, b: str) -> list[float]:
    for name in (a, b):
        if name not in rows[0]:
            raise SystemExit(f"Column {name!r} not found in the CSV")
    out = []
    for r in rows:
        ta, tb = (r[a] or "").strip(), (r[b] or "").strip()
        if ta and tb:
            out.append((float(ta) - float(tb)) * MM_TO_CM)
    return out


def summarise(vals: list[float]) -> dict:
    cuts = statistics.quantiles(vals, n=40, method="inclusive")  # 39 cut points, 2.5% steps
    return {"n": len(vals), "mean": statistics.fmean(vals), "sd": statistics.stdev(vals),
            "p2_5": cuts[0], "p97_5": cuts[38]}


def stats_for(rows, spec) -> dict | None:
    if spec is None:
        return None
    if spec[0] == "col":
        return summarise(column(rows, spec[1], spec[2]))
    return summarise(diff(rows, spec[1], spec[2]))


def q(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', "'") + '"'


def population(rows) -> dict:
    ages = column(rows, "Age", 1.0)
    dates = [r["Date"] for r in rows if r.get("Date")]
    return {"n": len(rows), "age_mean": statistics.fmean(ages), "age_min": min(ages), "age_max": max(ages),
            "first_date": min(dates, key=_date_key), "last_date": max(dates, key=_date_key)}


def _date_key(text: str):
    return datetime.strptime(text.strip(), "%d-%b-%y")  # e.g. 4-Oct-10


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    male_rows, female_rows = read(sys.argv[1]), read(sys.argv[2])
    pop_m, pop_f = population(male_rows), population(female_rows)

    lines = [
        "group: anthropometry",
        "description: >",
        "  Body size, lengths, breadths and girths for the adult male, computed from the public ANSUR II",
        f"  male file (n={pop_m['n']}, US Army personnel, mean age {pop_m['age_mean']:.1f}). This is a",
        "  military sample, not a general-population reference. Parameters without a direct ANSUR II",
        "  measurement have no value. Generated by tools/build_ansur_data.py; edit that script, not this file.",
        "parameters:",
    ]
    female_out = []
    for pid, name, unit, spec, note in P:
        m = stats_for(male_rows, spec)
        if m is None:
            extra = f", note: {q(note)}" if note else ""
            lines.append(f"  - {{id: {pid}, name: {q(name)}, unit: {unit if '/' not in unit else q(unit)}, kind: anchor{extra}}}")
            continue
        f = stats_for(female_rows, spec)
        n_txt = f"n={m['n']}"
        full_note = (note + " " if note else "") + f"ANSUR II male, {n_txt}."
        ukey = unit if "/" not in unit else q(unit)
        lines.append(
            f"  - {{id: {pid}, name: {q(name)}, unit: {ukey}, kind: anchor, "
            f"values: {{male: {{adult: {m['mean']:.2f}}}}}, sd: {{male: {{adult: {m['sd']:.2f}}}}}, "
            f"status: referenced, source_ref: ansur2, "
            f"ref_range: {{male: [{m['p2_5']:.2f}, {m['p97_5']:.2f}], female: [{f['p2_5']:.2f}, {f['p97_5']:.2f}]}}, "
            f"range_source_ref: ansur2, note: {q(full_note)}}}")
        female_out.append((pid, unit, f, m))

    OUT_VALUES.parent.mkdir(parents=True, exist_ok=True)
    OUT_VALUES.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fl = [
        "# Female ANSUR II statistics. These are the VALIDATION TARGETS for the Phase 2 sex rules:",
        "# the female derived from the male baseline must reproduce these within tolerance.",
        "# They are deliberately not loaded as baseline values. Generated by tools/build_ansur_data.py.",
        "sex: female",
        "source_ref: ansur2",
        f"population: {q('ANSUR II female, US Army personnel, n=%d, mean age %.1f' % (pop_f['n'], pop_f['age_mean']))}",
        "parameters:",
    ]
    for pid, unit, f, m in female_out:
        fl.append(
            f"  {pid}: {{unit: {unit if '/' not in unit else q(unit)}, mean: {f['mean']:.2f}, sd: {f['sd']:.2f}, "
            f"p2_5: {f['p2_5']:.2f}, p97_5: {f['p97_5']:.2f}, n: {f['n']}, ratio_to_male: {f['mean'] / m['mean']:.4f}}}")
    OUT_FEMALE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FEMALE.write_text("\n".join(fl) + "\n", encoding="utf-8")

    print(f"male:   {pop_m}")
    print(f"female: {pop_f}")
    print(f"wrote {len(female_out)} parameters with values, "
          f"{len(P) - len(female_out)} without an ANSUR II equivalent")


if __name__ == "__main__":
    main()
