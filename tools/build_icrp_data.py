"""Build the ICRP Publication 89 reference data files.

Usage:  python tools/build_icrp_data.py        (run from the repository root)

Writes:
    humansim/data/values/icrp_composition.yaml   organ and tissue masses by age and sex (Table 2.8)
    humansim/data/values/icrp_physiology.yaml    everything else (Tables 2.9, 2.11-2.14, 2.24, 2.27,
                                                 2.29-2.31, 2.39-2.43, section 7.4, ...)

The numbers below were transcribed from the ICRP 89 PDF supplied by the user (text layer of the
summary tables in chapter 2). The script refuses to write unless the transcription passes the
independent checks that the report itself provides:
  * the non-duplicated tissue rows of Table 2.8 sum to ~96% of total body mass (footnote c)
  * regional blood flow rates (Table 2.40, liver arterial only) sum to 100% of cardiac output
  * regional blood volumes (Table 2.14) sum to 100% of total blood volume
  * skeletal components of Table 2.8 add up to the total skeleton
Ages: 0 = newborn, 1, 5, 10, 15 years, "adult". Where the report gives one value for both sexes
below 15 years, it is stored for both sexes.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_C = ROOT / "humansim" / "data" / "values" / "icrp_composition.yaml"
OUT_P = ROOT / "humansim" / "data" / "values" / "icrp_physiology.yaml"
AGES = [0, 1, 5, 10, 15, "adult"]
SRC = "icrp89"
N = None  # not given in the report


def q(t: str) -> str:
    return '"' + t.replace("\\", "\\\\").replace('"', "'") + '"'


def fmt(x: float) -> str:
    return repr(round(x, 10))


def split(cell):
    """'1310/1180' -> (1310.0, 1180.0); 5 -> (5.0, 5.0); None -> (None, None)."""
    if cell is None:
        return None, None
    if isinstance(cell, str):
        a, b = cell.split("/")
        return float(a), float(b)
    return float(cell), float(cell)


def series(row8):
    """Table 2.8 layout [newborn, 1, 5, 10, 15 M, 15 F, adult M, adult F] -> {sex: {age: value}}."""
    nb, y1, y5, y10, m15, f15, ma, fa = row8
    out = {"male": {}, "female": {}}
    for age, cell in zip((0, 1, 5, 10), (nb, y1, y5, y10)):
        m, f = split(cell)
        out["male"][age], out["female"][age] = m, f
    out["male"][15], out["female"][15] = split(m15)[0], split(f15)[0]
    out["male"]["adult"], out["female"]["adult"] = split(ma)[0], split(fa)[0]
    return {s: {a: v for a, v in d.items() if v is not None} for s, d in out.items()}


def six(male6, female6=None):
    """Six values (newborn..adult) per sex; female defaults to the same as male."""
    female6 = female6 or male6
    return {"male": {a: v for a, v in zip(AGES, male6) if v is not None},
            "female": {a: v for a, v in zip(AGES, female6) if v is not None}}


def adult(m, f):
    return {"male": {"adult": m} if m is not None else {}, "female": {"adult": f} if f is not None else {}}


def scale(vals, k):
    return {s: {a: v * k for a, v in d.items()} for s, d in vals.items()}


# --- Table 2.8: organ and tissue masses (g). (id, name, dup_flag, row) ---------------------------
T28 = [
    ("adipose_tissue_mass", "Adipose tissue (total)", True, [930, 3800, 5500, 8600, 12000, 18700, 18200, 22500]),
    ("separable_adipose_mass", "Separable adipose tissue, excluding yellow marrow", False, [890, 3600, 5000, 7500, 9500, 16000, 14500, 19000]),
    ("adrenals_mass", "Adrenal glands (2)", False, [6, 4, 5, 7, 10, 9, 14, 13]),
    ("tongue_mass", "Tongue", False, [3.5, 10, 19, 32, 56, 53, 73, 60]),
    ("salivary_glands_mass", "Salivary glands", False, [6, 24, 34, 44, 68, 65, 85, 70]),
    ("oesophagus_wall_mass", "Oesophagus wall", False, [2, 5, 10, 18, 30, 30, 40, 35]),
    ("stomach_wall_mass", "Stomach wall", False, [7, 20, 50, 85, 120, 120, 150, 140]),
    ("stomach_contents_mass", "Stomach contents", False, [40, 67, 83, 117, 200, 200, 250, 230]),
    ("small_intestine_wall_mass", "Small intestine wall", False, [30, 85, 220, 370, 520, 520, 650, 600]),
    ("small_intestine_contents_mass", "Small intestine contents", False, [56, 93, 117, 163, 280, 280, 350, 280]),
    ("right_colon_wall_mass", "Right colon wall", False, [7, 20, 49, 85, 122, 122, 150, 145]),
    ("right_colon_contents_mass", "Right colon contents", False, [24, 40, 50, 70, 120, 120, 150, 160]),
    ("left_colon_wall_mass", "Left colon wall", False, [7, 20, 49, 85, 122, 122, 150, 145]),
    ("left_colon_contents_mass", "Left colon contents", False, [12, 20, 25, 35, 60, 60, 75, 80]),
    ("rectosigmoid_wall_mass", "Rectosigmoid wall", False, [3, 10, 22, 40, 56, 56, 70, 70]),
    ("rectosigmoid_contents_mass", "Rectosigmoid contents", False, [12, 20, 25, 35, 60, 60, 75, 80]),
    ("liver_mass", "Liver", False, [130, 330, 570, 830, 1300, 1300, 1800, 1400]),
    ("gallbladder_wall_mass", "Gallbladder wall", False, [0.5, 1.4, 2.6, 4.4, 7.7, 7.3, 10, 8]),
    ("gallbladder_contents_mass", "Gallbladder contents", False, [2.8, 8, 15, 26, 45, 42, 58, 48]),
    ("pancreas_mass", "Pancreas", False, [6, 20, 35, 60, 110, 100, 140, 120]),
    ("brain_mass", "Brain", False, [380, 950, "1310/1180", "1400/1220", 1420, 1300, 1450, 1300]),
    ("breasts_mass", "Breasts", False, [N, N, N, N, 15, 250, 25, 500]),
    ("heart_mass_with_blood", "Heart, with blood", True, [46, 98, 220, 370, 660, 540, 840, 620]),
    ("heart_mass", "Heart, tissue only", False, [20, 50, 85, 140, 230, 220, 330, 250]),
    ("blood_mass", "Blood", False, [290, 530, 1500, 2500, 4800, 3500, 5600, 4100]),
    ("eyes_mass", "Eyes (2)", False, [6, 7, 11, 12, 13, 13, 15, 15]),
    ("storage_fat_mass", "Fat (storage fat)", True, [370, 2300, 3600, 6000, 9000, 14000, 14600, 18000]),
    ("skin_mass", "Skin", False, [175, 350, 570, 820, 2000, 1700, 3300, 2300]),
    ("skeletal_muscle_mass", "Skeletal muscle", False, [800, 1900, 5600, 11000, 24000, 17000, 29000, 17500]),
    ("pituitary_mass", "Pituitary gland", False, [0.1, 0.15, 0.25, 0.35, 0.5, 0.5, 0.6, 0.6]),
    ("larynx_mass", "Larynx", False, [1.3, 4, 7, 12, 22, 15, 28, 19]),
    ("trachea_mass", "Trachea", False, [0.5, 1.5, 2.5, 4.5, 7.5, 6, 10, 8]),
    ("lungs_mass_with_blood", "Lungs, with blood", True, [60, 150, 300, 500, 900, 750, 1200, 950]),
    ("lungs_mass", "Lungs, tissue only", False, [30, 80, 125, 210, 330, 290, 500, 420]),
    ("skeleton_mass", "Skeleton (total)", True, [370, 1170, 2430, 4500, 7950, 7180, 10500, 7800]),
    ("cortical_bone_mass", "Bone, cortical", False, [135, 470, 1010, 1840, 3240, 2960, 4400, 3200]),
    ("trabecular_bone_mass", "Bone, trabecular", False, [35, 120, 250, 460, 810, 740, 1100, 800]),
    ("bone_mass_total", "Bone, total", True, [170, 590, 1260, 2300, 4050, 3700, 5500, 4000]),
    ("active_marrow_mass", "Marrow, active", False, [50, 150, 340, 630, 1080, 1000, 1170, 900]),
    ("inactive_marrow_mass", "Marrow, inactive", False, [0, 20, 160, 630, 1480, 1380, 2480, 1800]),
    ("cartilage_mass", "Cartilage", False, [130, 360, 600, 820, 1140, 920, 1100, 900]),
    ("teeth_mass", "Teeth", False, [0.7, 5, 15, 30, 45, 35, 50, 40]),
    ("skeleton_miscellaneous_mass", "Skeleton, miscellaneous", False, [20, 45, 55, 90, 155, 145, 200, 160]),
    ("spleen_mass", "Spleen", False, [9.5, 29, 50, 80, 130, 130, 150, 130]),
    ("thymus_mass", "Thymus", False, [13, 30, 30, "40/35", 35, 30, 25, 20]),
    ("thyroid_mass", "Thyroid", False, [1.3, 1.8, 3.4, 7.9, 12, 12, 20, 17]),
    ("tonsils_mass", "Tonsils (2 palatine)", False, [0.1, 0.5, 2, 3, 3, 3, 3, 3]),
    ("kidneys_mass", "Kidneys (2)", False, [25, 70, 110, 180, 250, 240, 310, 275]),
    ("ureters_mass", "Ureters (2)", False, [0.77, 2.2, 4.2, 7.0, 12, 12, 16, 15]),
    ("urinary_bladder_mass", "Urinary bladder", False, [4, 9, 16, 25, 40, 35, 50, 40]),
    ("urethra_mass", "Urethra", False, ["0.48/0.14", "1.4/0.42", "2.6/0.78", "4.4/1.3", 7.7, 2.3, 10, 3]),
    ("testes_mass", "Testes (2)", False, [0.85, 1.5, 1.7, 2, 16, N, 35, N]),
    ("epididymes_mass", "Epididymes (2)", False, [0.25, 0.35, 0.45, 0.60, 1.6, N, 4, N]),
    ("prostate_mass", "Prostate", False, [0.8, 1.0, 1.2, 1.6, 4.3, N, 17, N]),
    ("ovaries_mass", "Ovaries (2)", False, [0.3, 0.8, 2.0, 3.5, N, 6, N, 11]),
    ("fallopian_tubes_mass", "Fallopian tubes (2)", False, [0.25, 0.25, 0.35, 0.50, N, 1.1, N, 2.1]),
    ("uterus_mass", "Uterus", False, [4.0, 1.5, 3, 4, N, 30, N, 80]),
]
BODY_MASS = [3.5, 10, 19, 32, 56, 53, 73, 60]  # Table 2.8 last row / Table 2.9, kg
HEIGHT = [51, 76, 109, 138, 167, 161, 176, 163]
BSA = [0.24, 0.48, 0.78, 1.12, 1.62, 1.55, 1.90, 1.66]
# Sex-specific organs: below the age where the report starts giving them for one sex, the other sex
# has no entry (testes for females and so on). The single-column early values of these organs apply
# to the sex that has the organ only.
MALE_ONLY = {"testes_mass", "epididymes_mass", "prostate_mass"}
FEMALE_ONLY = {"ovaries_mass", "fallopian_tubes_mass", "uterus_mass"}

# --- regional flow and blood content (Tables 2.40, 2.14), % ------------------------------------------
# (id suffix, name, male %, female %)
FLOW = [
    ("fat", "Fat", 5.0, 8.5), ("brain", "Brain", 12, 12), ("stomach_oesophagus", "Stomach and oesophagus", 1.0, 1.0),
    ("small_intestine", "Small intestine", 10, 11), ("large_intestine", "Large intestine", 4.0, 5.0),
    ("coronary", "Coronary tissue", 4.0, 5.0), ("kidneys", "Kidneys", 19, 17),
    ("bronchial", "Bronchial tissue", 2.5, 2.5), ("skeletal_muscle", "Skeletal muscle", 17, 12),
    ("pancreas", "Pancreas", 1.0, 1.0), ("skeleton", "Skeleton (total)", 5.0, 5.0),
    ("skin", "Skin", 5.0, 5.0), ("spleen", "Spleen", 3.0, 3.0), ("thyroid", "Thyroid", 1.5, 1.5),
    ("lymph_nodes", "Lymph nodes", 1.7, 1.7), ("gonads", "Gonads", 0.05, 0.02),
    ("adrenals", "Adrenals", 0.3, 0.3), ("urinary_bladder", "Urinary bladder", 0.06, 0.06),
    ("other_tissues", "All other tissues", 1.39, 1.92),
]
FLOW_SKELETON_PARTS = [("red_marrow", "Red marrow", 3.0, 3.0), ("trabecular_bone", "Trabecular bone", 0.9, 0.9),
                       ("cortical_bone", "Cortical bone", 0.6, 0.6), ("other_skeleton", "Other skeleton", 0.5, 0.5)]
LIVER_ARTERIAL = (6.5, 6.5)
LIVER_TOTAL = (25.5, 27.0)
BLOOD = [
    ("fat", "Fat", 5.0, 8.5), ("brain", "Brain", 1.2, 1.2), ("stomach_oesophagus", "Stomach and oesophagus", 1.0, 1.0),
    ("small_intestine", "Small intestine", 3.8, 3.8), ("large_intestine", "Large intestine", 2.2, 2.2),
    ("right_heart", "Right heart", 4.5, 4.5), ("left_heart", "Left heart", 4.5, 4.5),
    ("coronary", "Coronary tissue", 1.0, 1.0), ("kidneys", "Kidneys", 2.0, 2.0), ("liver", "Liver", 10, 10),
    ("pulmonary", "Pulmonary", 10.5, 10.5), ("bronchial", "Bronchial tissue", 2.0, 2.0),
    ("skeletal_muscle", "Skeletal muscle", 14, 10.5), ("pancreas", "Pancreas", 0.6, 0.6),
    ("skeleton", "Skeleton (total)", 7.0, 7.0), ("skin", "Skin", 3.0, 3.0), ("spleen", "Spleen", 1.4, 1.4),
    ("thyroid", "Thyroid", 0.06, 0.06), ("lymph_nodes", "Lymph nodes", 0.2, 0.2), ("gonads", "Gonads", 0.04, 0.02),
    ("adrenals", "Adrenals", 0.06, 0.06), ("urinary_bladder", "Urinary bladder", 0.02, 0.02),
    ("other_tissues", "All other tissues", 1.92, 1.92), ("aorta_large_arteries", "Aorta and large arteries", 6.0, 6.0),
    ("large_veins", "Large veins", 18, 18),
]


def check_transcription() -> list[str]:
    msgs = []
    # 1. Table 2.8 sum to ~96% of body mass, per column
    cols = ["newborn", "1 y", "5 y", "10 y", "15 y M", "15 y F", "adult M", "adult F"]
    for i, col in enumerate(cols):
        total = 0.0
        for pid, _, dup, row in T28:
            if dup:
                continue
            cell = row[i]
            if cell is None:
                continue
            total += split(cell)[0 if i not in (0, 1, 2, 3) else 0] if i < 4 and False else _cell_value(cell, i, pid)
        frac = total / 1000 / BODY_MASS[i]
        msgs.append((col, total / 1000, BODY_MASS[i], frac))
        if not 0.93 <= frac <= 0.99:
            raise SystemExit(f"Table 2.8 column {col}: listed rows sum to {frac:.1%} of body mass (expected ~96%)")
    # 2. skeleton parts add up to total skeleton
    ix = {pid: row for pid, _, _, row in T28}
    for i, col in enumerate(cols):
        parts = ["cortical_bone_mass", "trabecular_bone_mass", "active_marrow_mass", "inactive_marrow_mass",
                 "cartilage_mass", "teeth_mass", "skeleton_miscellaneous_mass"]
        s = sum(split(ix[p][i])[0] for p in parts)
        t = split(ix["skeleton_mass"][i])[0]
        if abs(s - t) / t > 0.03:
            raise SystemExit(f"skeleton parts {s} != total {t} in column {col}")
    # 3. flows and blood content
    for k, name in ((2, "male"), (3, "female")):
        flow = sum(r[k] for r in FLOW) + LIVER_ARTERIAL[k - 2]
        blood = sum(r[k] for r in BLOOD)
        if abs(flow - 100) > 0.1:
            raise SystemExit(f"regional flow ({name}) sums to {flow}, not 100")
        if abs(blood - 100) > 0.05:
            raise SystemExit(f"regional blood content ({name}) sums to {blood}, not 100")
    return msgs


def _cell_value(cell, col_index, pid):
    a, b = split(cell)
    # columns: 0-3 same for both sexes except M/F pairs (take male); 4 = 15 M; 5 = 15 F; 6 = adult M; 7 = adult F
    return a


# --- YAML emitters ----------------------------------------------------------------------------------
def values_yaml(vals, scale_k=1.0):
    parts = []
    for sex in ("male", "female"):
        d = vals.get(sex) or {}
        if not d:
            continue
        inner = ", ".join(f"{a}: {fmt(v * scale_k)}" for a, v in d.items())
        parts.append(f"{sex}: {{{inner}}}")
    return "{" + ", ".join(parts) + "}"


def record(pid, name, unit, vals, note, kind="anchor", scale_k=1.0, source=SRC):
    if not any(vals.values()):
        raise SystemExit(f"{pid}: no values")
    return (f"  - {{id: {pid}, name: {q(name)}, unit: {q(unit) if '/' in unit or '^' in unit else unit}, kind: {kind}, "
            f"values: {values_yaml(vals, scale_k)}, status: referenced, source_ref: {source}, note: {q(note)}}}")


def build():
    comp, phys = [], []
    sex_note = "ICRP 89 reference man (73 kg) and woman (60 kg); one value for both sexes below 15 years where the report gives one."
    dup_note = " The report marks this entry as duplicating other rows, so it is excluded from whole-body sums."
    for pid, name, dup, row in T28:
        vals = series(row)
        if pid in MALE_ONLY:
            vals["female"] = {}
            vals["male"] = {a: v for a, v in vals["male"].items()}
        if pid in FEMALE_ONLY:
            vals["male"] = {}
        comp.append(record(pid, name + " mass" if not name.lower().endswith("mass") and not pid.endswith("_mass") else name,
                           "kg", vals, "ICRP 89 Table 2.8. " + sex_note + (dup_note if dup else ""), scale_k=1e-3))

    cb = [BODY_MASS[0:4] + [BODY_MASS[4]], None]
    ref = series(BODY_MASS)
    phys.append(record("reference_body_mass", "Reference person body mass (ICRP 89)", "kg", ref,
                       "ICRP 89 Tables 2.8 and 2.9. The mass of the person the organ and physiology data describe; differs from the ANSUR II body mass."))
    phys.append(record("reference_height", "Reference person height (ICRP 89)", "cm", series(HEIGHT), "ICRP 89 Table 2.9."))
    phys.append(record("reference_body_surface_area", "Reference person body surface area (ICRP 89)", "m^2", series(BSA), "ICRP 89 Tables 2.9 and 2.26."))

    phys.append(record("blood_volume", "Total blood volume", "L", six([0.27, 0.5, 1.4, 2.4, 4.5, 5.3], [0.27, 0.5, 1.4, 2.4, 3.3, 3.9]),
                       "ICRP 89 section 7.4, reference values for total blood volume.", kind="regulated"))
    phys.append(record("cardiac_output", "Cardiac output (resting)", "L/min", six([0.6, 1.2, 3.4, 5.0, 6.1, 6.5], [0.6, 1.2, 3.4, 5.0, 6.1, 5.9]),
                       "ICRP 89 Table 2.39: trimmed means from >50 studies of healthy resting recumbent subjects.", kind="regulated"))
    bmr_m, bmr_f = [8.5, 24, 45, 53, 63, 68], [8.5, 23, 43, 48, 51, 52]
    phys.append(record("resting_metabolic_rate", "Basal metabolic rate", "kcal/day",
                       six([x * 24 for x in bmr_m], [x * 24 for x in bmr_f]),
                       "ICRP 89 Table 2.29, basal metabolic rate in kcal/h multiplied by 24.", kind="regulated"))
    phys.append(record("total_energy_expenditure", "Total energy expenditure", "kcal/day",
                       six([500, 800, 1600, 1900, 2400, 2800], [500, 750, 1400, 1700, 1800, 1800]), "ICRP 89 Table 2.29."))
    lf = LIVER_TOTAL
    phys.append(record("hepatic_flow_fraction", "Hepatic blood flow, fraction of cardiac output", "dimensionless",
                       adult(lf[0] / 100, lf[1] / 100), "ICRP 89 Table 2.40, liver total (arterial 6.5% plus portal).", kind="regulated"))
    phys.append(record("hepatic_arterial_flow_fraction", "Hepatic arterial flow, fraction of cardiac output", "dimensionless",
                       adult(LIVER_ARTERIAL[0] / 100, LIVER_ARTERIAL[1] / 100), "ICRP 89 Table 2.40."))
    for suf, name, m, f in FLOW + FLOW_SKELETON_PARTS:
        phys.append(record(f"blood_flow_fraction_{suf}", f"Blood flow to {name.lower()}, fraction of cardiac output", "dimensionless",
                           adult(m / 100, f / 100), "ICRP 89 Table 2.40. Adult reference values."))
    for suf, name, m, f in BLOOD:
        phys.append(record(f"blood_content_fraction_{suf}", f"Blood content of {name.lower()}, fraction of total blood volume", "dimensionless",
                           adult(m / 100, f / 100), "ICRP 89 Table 2.14. Adult reference values."))
    resp = [("total_lung_capacity", "Total lung capacity", [0.28, 0.55, 1.6, 2.9], [5.4, 7.0], [4.5, 5.0]),
            ("functional_residual_capacity", "Functional residual capacity", [0.15, 0.24, 0.77, 1.5], [2.7, 3.3], [2.3, 2.7]),
            ("vital_capacity", "Vital capacity", [0.20, 0.38, 1.0, 2.3], [4.0, 5.0], [3.3, 3.6]),
            ("anatomical_dead_space", "Anatomical dead space", [0.014, 0.020, 0.046, 0.078], [0.13, 0.15], [0.11, 0.12])]
    for pid, name, kid, m, f in resp:
        ages = [0.25, 1, 5, 10, 15, "adult"]
        phys.append(record(pid, name, "L", {"male": dict(zip(ages, kid + m)), "female": dict(zip(ages, kid + f))},
                           "ICRP 89 Table 2.31 (rounded values from ICRP Publication 66). Ages: 3 months, 1, 5, 10, 15 years, adult."))
    phys.append(record("water_requirement", "Water intake in food and fluids", "L/day", adult(2.6, 1.96), "ICRP 89 Table 2.30.", kind="regulated"))
    phys.append(record("water_from_food_oxidation", "Water from oxidation of food", "L/day", adult(0.3, 0.225), "ICRP 89 Table 2.30."))
    phys.append(record("water_loss_sweat", "Water loss in sweat", "L/day", adult(0.5, 0.375), "ICRP 89 Table 2.30."))
    phys.append(record("water_loss_insensible", "Insensible water loss (lungs and skin)", "L/day", adult(0.69, 0.515), "ICRP 89 Table 2.30."))
    phys.append(record("water_loss_faeces", "Water loss in faeces", "L/day", adult(0.11, 0.095), "ICRP 89 Table 2.30."))
    phys.append(record("urinary_excretion_volume", "Urine volume", "mL/day", six([300, 400, 500, 700, 1200, 1600], [300, 400, 500, 700, 1200, 1200]),
                       "ICRP 89 Table 2.41."))
    phys.append(record("creatinine_excretion", "Urinary creatinine excretion", "g/day", six([0.05, 0.11, 0.33, 0.65, 1.4, 1.7], [0.05, 0.11, 0.33, 0.65, 1.0, 1.0]),
                       "ICRP 89 Table 2.42."))
    phys.append(record("cortical_bone_density", "Cortical bone density (hydrated)", "g/cm^3", six([1.65, 1.66, 1.70, 1.75, 1.80, 1.90]),
                       "ICRP 89 Table 2.20, hydrated cortical bone."))
    phys.append(record("cortical_bone_remodelling_rate", "Cortical bone remodelling rate", "1/year", six([3.0, 1.05, 0.56, 0.33, 0.19, 0.03]).__class__({
        "male": dict(zip([0.25, 1, 5, 10, 15, "adult"], [3.0, 1.05, 0.56, 0.33, 0.19, 0.03])),
        "female": dict(zip([0.25, 1, 5, 10, 15, "adult"], [3.0, 1.05, 0.56, 0.33, 0.19, 0.03]))}),
        "ICRP 89 Table 2.43 (300, 105, 56, 33, 19, 3 percent per year). Ages: 0-3 months, 1, 5, 10, 15 years, adult."))
    phys.append(record("trabecular_bone_remodelling_rate", "Trabecular bone remodelling rate", "1/year", {
        "male": dict(zip([0.25, 1, 5, 10, 15, "adult"], [3.0, 1.05, 0.66, 0.48, 0.35, 0.18])),
        "female": dict(zip([0.25, 1, 5, 10, 15, "adult"], [3.0, 1.05, 0.66, 0.48, 0.35, 0.18]))},
        "ICRP 89 Table 2.43 (300, 105, 66, 48, 35, 18 percent per year). Ages: 0-3 months, 1, 5, 10, 15 years, adult."))
    phys.append(record("skeletal_calcium_mass", "Skeletal calcium", "kg",
                       series([28, 100, 240, 460, 830, 760, 1180, 860]), "ICRP 89 Table 2.24.", scale_k=1e-3))
    phys.append(record("epidermis_mass", "Epidermis mass", "kg", series([12, 24, 39, 56, 100, 80, 120, 85]), "ICRP 89 Table 2.27.", scale_k=1e-3))
    phys.append(record("dermis_mass", "Dermis mass", "kg", series([163, 326, 531, 764, 1900, 1620, 3180, 2215]), "ICRP 89 Table 2.27.", scale_k=1e-3))
    phys.append(record("lymphocyte_mass", "Total lymphocyte mass", "kg", six([0.15, 0.65, 0.65, 0.9, 1.25, 1.5], [0.15, 0.65, 0.65, 0.9, 1.25, 1.3]),
                       "ICRP 89 Table 2.16."))
    phys.append(record("fixed_lymphatic_tissue_mass", "Fixed lymphatic tissue mass", "kg", adult(0.73, 0.60), "ICRP 89 Table 2.15."))
    seg = [("oesophagus_length", "Oesophagus length", [10, 13, 18, 23, 27, 28], [10, 13, 18, 23, 26, 26]),
           ("small_intestine_length", "Small intestine length", [80, 120, 170, 220, 270, 280], [80, 120, 170, 220, 260, 260]),
           ("right_colon_length", "Right colon length", [14, 18, 23, 28, 35, 34], [14, 18, 23, 28, 30, 30]),
           ("left_colon_length", "Left colon length", [16, 21, 26, 31, 35, 38], [16, 21, 26, 31, 35, 35]),
           ("rectosigmoid_length", "Rectosigmoid length", [15, 21, 26, 31, 35, 38], [15, 21, 26, 31, 35, 35]),
           ("large_intestine_length", "Large intestine length (total)", [45, 60, 75, 90, 100, 110], [45, 60, 75, 90, 100, 100])]
    for pid, name, m, f in seg:
        phys.append(record(pid, name, "cm", six(m, f), "ICRP 89 Table 2.11."))
    phys.append(record("bone_tissue_volume", "Volume of bone tissue", "mL", adult(2710, None), "ICRP 89 Table 2.19 (adult male only)."))
    phys.append(record("cortical_bone_volume", "Volume of cortical bone", "mL", adult(2130, None), "ICRP 89 Table 2.19 (adult male only)."))
    phys.append(record("trabecular_bone_volume", "Volume of trabecular bone", "mL", adult(580, None), "ICRP 89 Table 2.19 (adult male only)."))
    return comp, phys


def main() -> None:
    sums = check_transcription()
    comp, phys = build()
    OUT_C.write_text(
        "group: composition\ndescription: >\n  Organ and tissue masses by age and sex from ICRP Publication 89 Table 2.8 (reference man 73 kg,\n"
        "  reference woman 60 kg). Generated by tools/build_icrp_data.py; edit that script, not this file.\n"
        "parameters:\n" + "\n".join(comp) + "\n", encoding="utf-8")
    OUT_P.write_text(
        "group: physiology\ndescription: >\n  Physiological and body reference values by age and sex from ICRP Publication 89 (chapter 2 summary\n"
        "  tables). Generated by tools/build_icrp_data.py; edit that script, not this file.\n"
        "parameters:\n" + "\n".join(phys) + "\n", encoding="utf-8")
    for col, tot, bm, frac in sums:
        print(f"Table 2.8 {col:>8}: listed tissues {tot:6.2f} kg of {bm} kg = {frac:.1%}")
    print(f"wrote {len(comp)} composition and {len(phys)} physiology parameters")


if __name__ == "__main__":
    main()
