# HumanSim

A Python library that models the human body as a network of interconnected anatomical and
physiological parameters, for population averages (not individuals). You set a context, edit any
value, and see what changes. See `docs/SPEC.md` for the full design.

**Status: all eight phases have a first version (core, sex, solver, development, hormones, pharmacology, genetics and gene edits, UI).** The core, derived layer, consistency checks, modifier machinery, adult male and
female baselines (both stored, either derivable from the other), a steady-state solver and an age
axis for birth to 20 years, sex hormones (levels, one dose-response effect, two dynamic models) exist. A first pharmacokinetics layer and a one-gene pharmacogenetics layer exist. A local web UI is included. Run `python -m humansim coverage` to see exactly how much data there is.

## Install and run

```
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest

python -m humansim show                                  # resolved adult male
python -m humansim show --set body_mass=100 --changed    # what changes if mass is 100 kg
python -m humansim show --set body_mass="220 lb" --changed
python -m humansim coverage --missing                    # what still needs a value
python -m humansim check                                 # consistency checks

python -m humansim show --sex female --trace             # derived female, with every rule applied
python -m humansim compare                               # male and female side by side
python -m humansim validate                              # rule/pair consistency and an independent NCHS check

python -m humansim show --set liver_mass=9 --solve --changed --trace   # the ripple, step by step
python -m humansim show --set height="7 ft" --solve --changed
python -m humansim check --solve --set liver_mass=9

python -m humansim show --age 10 --sex female            # a 10-year-old girl (CDC growth curves)
python -m humansim show --age 0.5 --z 1                  # 6 months, one SD above the median
python -m humansim show --age 3 --solve --set height=100
```

```python
from humansim import Registry, Context, resolve

reg = Registry.load()
state = resolve(reg, Context(edits={"body_mass": 100.0}))
state.get("body_mass_index")      # derived, recomputed
state.flags["body_mass"]          # {'edited'}
reg["height"].spread("male", "adult")   # standard deviation, where the data has one
```

## How values are marked

Every stored parameter has a status:

| Status | Meaning |
|---|---|
| `referenced` | Value computed from, or read in, a source named in `data/sources/sources.yaml` |
| `unsourced` | A value is present but I could not check it against a source |
| `placeholder` | A stand-in |
| `missing` | The parameter is defined but has no value yet |

`derived` values are computed from stored values and are never stored. They inherit `unsourced` and
`edited` flags from their inputs. Where a source gives one, a parameter also carries a standard
deviation and a 2.5th to 97.5th percentile reference range; values edited outside that range are
flagged `outside_reference_range`.

## What is real and what is not

- **ANSUR II (computed from the raw files, 32 parameters, both sexes):** height, body mass, and the
  lengths, breadths and girths that ANSUR II measures, with SD and 2.5-97.5 percentile range.
  `tools/build_ansur_data.py` regenerates `data/values/anthropometry.yaml` from the public CSVs.
  **Caveat: the sample is US Army personnel (n=4,082 men, mean age 30), not the general public.**
  Their mean mass is 85.5 kg (men) and 67.8 kg (women).
- **ICRP Publication 89 (read from the PDF, 142 parameters, both sexes, birth to adult):** organ and
  tissue masses (Table 2.8, 57 rows), reference height, mass and surface area, blood volume, cardiac
  output, basal metabolic rate, energy expenditure, regional blood flow and blood content, respiratory
  volumes, water balance, urine and creatinine, bone density and remodelling, gut lengths, lymphatic
  masses. `tools/build_icrp_data.py` holds the transcription and **refuses to write unless the report's
  own totals check out**: tissue rows sum to 95.5-96.7% of body mass in every column (the report says
  96%), regional flows sum to 100% of cardiac output, regional blood contents to 100%, skeletal parts to
  the total skeleton. The report's reference people are a 73 kg man and a 60 kg woman, mainly
  European and North American data: a different population from ANSUR II.
- **Corrected by the PDF:** female adipose tissue is 22.5 kg (I had recalled 18.0), female skin 2.3 kg
  (1.7), female hepatic blood flow 27% of cardiac output (I had assumed the male 25.5%), and the
  lungs value I used (1.2 kg) was lungs *with blood*; tissue only is 0.5 kg.
- **Secondary or unsourced (6):** haematocrit (secondary page; ICRP gives none), core temperature,
  haemoglobin, SaO2, PaO2, brain power, muscle density: typical figures not retrieved.
- **Missing (22):** no ANSUR II equivalent (femur and tibia length, knee girth, standing reach,
  head height, neck length, body density...) or no source yet (enamel, bite force...).
- **Two populations in one model.** Body dimensions are ANSUR II (US Army: 85.5 kg men, 67.8 kg
  women); organs and physiology are the ICRP reference people (73 kg, 60 kg). The organ budget,
  metabolic cross-checks and the Mifflin check therefore compare against the ICRP reference mass when
  the state is unedited, and against `body_mass` once something is edited or solved. The solver works
  in ratios against the stored baseline, so it does not hide this and does not change stored values.
- **Left out on purpose:** silk, rete, structural lipid and the cold and working energy figures in
  the original document, which describe a non-human body plan.

## Sex: both directions

Neither sex is derived from the other by default. Each sex has its own stored value wherever data
exists: ANSUR II means, SDs and ranges for both sexes (32 parameters), ICRP 89 female and male
organ masses, blood volume, cardiac output, metabolic rate, flows and more (the table values are given
per sex, so most rules are documentation).

Sex rules (`data/modifiers/sex_female.yaml`, generated by `tools/build_sex_rules.py` from the stored
pairs) record `female = male * ratio` and are used in **both directions, only to fill a gap**:

- only a male value stored: female = male * ratio
- only a female value stored: male = female / ratio (the trace says "inverse of female rule")
- both stored: the rule is documentation; stored values always win, so nothing is applied twice

A value with no rule for the missing sex is left empty and flagged `no_sex_rule` (today only
`brain_metabolic_load`); `--keep-unadjusted` keeps the other sex's value, still flagged. Rule types:
`scaling` (size alone explains it, within 3%), `anatomical` (shape difference), `physiological`
(measured, no cause attributed, or tied to another value such as haemoglobin to haematocrit),
`assumed_neutral` (no sex difference assumed, unsourced).

**Honest limits.** `humansim validate` first checks that rules agree with the stored pairs (true by
construction; it catches stale files). Its second part is the independent one: female/male ratios
from NCHS general-population data. Height agrees to 0.6%, waist to 3.5%, but body mass differs by
8%: the Army sample has a lower female/male mass ratio than the US public. Treat ANSUR II sex ratios
as military-sample ratios.

## The solver (Phase 3)

`solve(registry, context, modifiers, couplings)` (or `--solve`) lets edits ripple. Couplings are in
`data/couplings/steady_state.yaml`. The ideas:

- **Reference-calibrated.** Every coupling returns a ratio and the solved value is `reference value *
  ratio`, where the reference is the same context with no edits. The unedited model reproduces its
  stored values exactly (tested), and a coupling only says how a value *moves*.
- **Accounting vs compensation.** Accounting couplings are bookkeeping and apply at once: body mass
  includes an edited organ, metabolic rate is the sum of tissue rates (Wang 2010), dimensions follow
  height. Compensation couplings are the body adjusting and have a **deadband**: a change inside
  normal variation moves nothing; beyond it the response starts from the band edge (continuous).
  Deadbands are either a per-coupling fraction (assumed 10%) or the driver's own 2.5-97.5 percentile range.
- **Labelled.** Each coupling is `mechanistic` (Fick principle, metabolic sum), `empirical` or
  `assumed`. Values moved by an assumed coupling carry `coupling_unsourced`; a driver outside the
  coupling's assumed validity range adds `outside_coupling_validity` (extrapolation).
- **Pinned edits and feedback.** A value you edit is never overwritten. Feedback loops (bigger
  heart, more metabolic demand, more output, bigger heart) are iterated to a fixed point; if the
  loop gain is 1 or more the result says `converged: False`.

Chains implemented: liver mass to body mass, metabolic rate, cardiac output (Fick), heart mass and
hepatic flow share; haemoglobin to arterial O2 content to cardiac output and heart mass; body mass to
organ, tissue and blood volume (outside the reference range only); height to every linear dimension
(isometric, head at power 0.5) and body mass (power 2, BMI-preserving).

Example, liver 1.8 to 9 kg (male): body mass +7.2 kg, metabolic rate 1808 to 3274 kcal/day, cardiac
output 6.5 to 10.7 L/min, heart 0.33 to 0.49 kg, hepatic flow 1.7 to 6.8 L/min (flagged as
extrapolation). A cross-check (`rmr_organ_crosscheck`) compares the organ-sum metabolic rate with the
Mifflin-St Jeor value; they agree within 3% for the male baseline and about 6% for the female.

**Not built, and why.** Liver to abdominal volume to lung compression, heat and water: there are no
abdominal-volume, lung-capacity or water data in the registry, and the solver does not guess.
Organ and tissue exponents (0.75, 1.0), the 10% deadbands, the validity ranges and the height
exponents are modelling assumptions, not fitted or sourced. Which tissue takes up a body-mass change
(fat or muscle) depends on cause and is not modelled.

## Development (Phase 4)

`Context(age=...)` takes `"adult"` or years. Below 20 years the model uses the **CDC 2000 growth
curves** (LMS tables in `data/growth/`): height, body mass and (to 36 months) head circumference, by
sex, interpolated in months and flagged `growth_curve` and `interpolated`. `growth_z` picks a
percentile (0 is the median, "the average child"). Reference ranges at a pediatric age are the 3rd to
97th percentile at that age, replacing the adult ranges.

Rules that keep pediatric and adult data apart:
- Height, body mass and head circumference (to 36 months) come from the CDC curves. Everything in ICRP 89
  has anchor ages (birth, 1, 5, 10, 15 years, adult): an exact anchor age gives the stored value, in
  between is **linear and flagged `interpolated`**, and 15 years to adult interpolates towards age 20.
- A parameter with neither a curve nor anchors **has no value** for a child (flagged `missing`); adult-only
  ICRP values (regional flows, hepatic flow fraction, haematocrit) are never used for a child.
- Sex rules are not used at pediatric ages; the ICRP tables give both sexes directly.
- Age 20 and over is the adult baseline, with a note that ageing after 20 is not modelled.
- The solver runs at any age; in a child deadbands use the age-specific 3rd-97th percentile ranges.

Checks: the LMS formula reproduces the published P3/P50/P97 from L, M, S in every table; BMI computed
from median height and weight agrees with the CDC median BMI within 8%; `check` reports the gap
between the CDC median at 20 years and the adult baseline. That gap is large for body mass
(men -17.5%, women -14%) because ANSUR II is a heavier military sample, and it is not hidden.

**Limits.** Only height, mass and head size have continuous curves; ICRP organs and physiology have
five anchor ages, so values between them are straight-line guesses. The organ-sum metabolic check is
skipped for children because the tissue rates (Wang 2010) are adult rates and children's organs burn
more per kg (the check gives -26% at age 10 if forced). Puberty timing, ageing after 20 and the
fetal tables (ICRP Tables 2.1-2.7) are not modelled. The CDC and WHO sites were unreachable; the CDC
tables came from a third-party PyPI copy and were spot-checked, not compared row by row. CDC data are
US children measured mainly 1963-1994.

## Hormones and dynamics (Phase 5)

**Levels.** `data/values/hormones.yaml` holds testosterone, SHBG, LH, FSH, oestradiol, progesterone,
cortisol and prolactin. Women's values depend on cycle phase (`Context(cycle_phase=...)`:
follicular, midcycle, luteal, postmenopausal; default follicular). A value here is the **midpoint of a
published single-laboratory reference interval, not a population mean**, and is flagged as such.

**Setting a hormone.** `Context(hormone_state={"testosterone_total": 1345})` or
`--hormone testosterone_total=1345`. The level is recorded in the trace and flagged `hormone_set`.

**One sourced effect.** Testosterone changes skeletal muscle mass by `2.95 kg x ln(T / T_ref)`
(`log_slope` modifier; slope fitted to the three dose groups of Bhasin 2001, fat-free mass taken as
muscle). It applies to men only and is flagged outside 542-2370 ng/dL, the studied range. With `--solve`
the change is pinned (couplings do not undo it) and ripples to body mass, metabolic rate, BMI and the
composition fractions. The solver's reference state now excludes hormone levels, so a hormone effect
counts as a change to respond to.

**Dynamics.** `dynamics.py` is a generic RK4 engine. `hpg.py` is a three-variable male feedback model
(LH, testicular responsiveness, testosterone) that sits at the reference with no input, is suppressed
by exogenous testosterone and recovers over weeks. `cycle.py` interpolates the stored phase values into a
28-day profile (`--sex female --cycle-day 14`).

**Limits.** The HPG structure, feedback strength and all three time constants are **assumptions**
(results carry `assumed`), not fitted to data. The cycle profile joins phase midpoints linearly, so it
misses the sharp LH peak. The Bhasin effect uses only the slope, from men 18-35 at 20 weeks of
suppressed-own-production dosing, and says nothing about haemoglobin or other tissues. No female feedback
model, no pediatric or ageing hormone levels yet.

## Pharmacology (Phase 6)

`python -m humansim pk --drug midazolam --dose 5 --sex female --set liver_mass=2.5 --solve`

`pharmacology.py` computes the plasma concentration of a drug after intravenous doses (boluses or
infusions, repeated if needed) in a one-compartment model. Each drug in `data/drugs/` carries a
reference clearance and volume per kg from a cited label. For a person they are scaled by that
person's own model values (liver mass, hepatic blood flow, body mass) relative to the same sex and age
unedited, so editing or solving the body changes the curve. Results are flagged `assumed_scaling`,
`one_compartment` and `scaled_from_changed_body`. Time is in hours.

Seed drugs: **acetaminophen IV** (FDA label: CL 0.27 L/h/kg, Vss 0.8 L/kg) and **midazolam IV** (Pfizer
label: CL 0.25-0.54 L/h/kg and V 1.0-3.1 L/kg; the model uses range midpoints, which are not means).

**Limits.** This is not PBPK. A single compartment understates the early peak (an acetaminophen 1 g
infusion gives about 16 mg/L against the label's 28) and gives half-lives 14% short (acetaminophen) and 20% long
(midazolam) of the label; `check` reports both. The scaling exponents are assumptions. Not modelled: oral
absorption (no sourced bioavailability or absorption rate), renal drugs (no sourced clearance values; the
gentamicin and caffeine pages I could reach had none), metabolites, protein binding, drug effects on
the body, and between-person variability. Per-kg label values are applied to both sexes because the labels
do not give volunteer sex or mass.

## Genetics (Phase 7)

`python -m humansim genes` lists genes, activity scores and population frequencies.
`python -m humansim pk --drug ibuprofen --dose 800 --infusion 1 --genotype CYP2C9='*1/*3'`

`genetics.py` implements the CPIC activity-score scheme: each allele has an activity value, a
diplotype's score is the sum, and the score gives a phenotype (normal, intermediate, poor). Set a
genotype with `Context(genotype={"CYP2C9": "*1/*3"})`. Allele frequencies give Hardy-Weinberg diplotype
and phenotype frequencies (`genes` command). A drug can list a `genotype_effects` entry: a table of
mean AUC ratios by activity score, and the model divides clearance by the ratio. The reference person is
treated as a normal metabolizer (*1/*1).

Built: CYP2C9 (*1, *2, *3; CPIC 2020 activity values) and **IV ibuprofen** (AUC x1.43 for *1/*3, from CPIC's
meta-analysis; the label's 800 mg AUC of 192.8 is reproduced for *1/*1).

**Limits.** One gene, one drug with genotype data. Only activity scores 2 and 1 have ibuprofen data:
1.5 is interpolated and flagged, and poor metabolizers (score below 1) are **refused**, not guessed. The
ratios come from oral studies, so applying them to IV dosing assumes bioavailability does not depend on
genotype. Allele frequencies exist for people of European ancestry only (about 11% *2 and 7% *3 from the
warfarin label); other alleles are counted as *1. CPIC and the label pages were read through a
fetch tool's summary, not page by page; verify the numbers against the originals before relying on
them. Ibuprofen clearance and volume are *derived* from the label's AUC and rate constant, with an assumed
73 kg volunteer. Genotype changes no anatomical or physiological parameter: no sourced numbers for that are
in the project yet. CYP3A5 was examined for midazolam and left out because the sources found show no effect.

## Gene edits: somatic and germline (sandbox)

```
python -m humansim catalogue
python -m humansim show --edit EPOR:truncated_gain_of_function --edit-mode germline --solve --changed
python -m humansim show --edit MSTN:loss_of_function --edit-fraction 0.5 --edit-tissue "skeletal muscle" --solve --changed
python -m humansim show --hypothesis FOXP2:heart_mass=ratio:1.3
python -m humansim pk --drug acetaminophen --dose 1000 --infusion 0.25 --edit MSTN:loss_of_function --solve
```

**This is a place to explore ideas, not a predictor.** Every edit effect has a tier. *Measured* means a
sourced effect size with an evidence grade. *Hypothetical* means a stated guess, or your own idea typed
in with `--hypothesis`; it is flagged `hypothetical` on every value it changes, on derived values built
from them, and as `downstream_of_hypothetical` on whatever the solver then moves. The CLI prints a
banner whenever a state contains one. Pharmacology results scaled from such a body carry
`hypothetical_inputs`.

- **Somatic** (in vivo): `--edit-fraction` is the share of cells edited. The effect is scaled linearly by it,
  which is an assumption (flag `mosaic_assumed_linear`). Edited tissue is recorded, not modelled.
- **Germline**: every cell, from conception (fraction 1, flag `germline`). Heritable, and prohibited or
  restricted for human embryos in many countries. Effects are applied at any age as a ratio on that
  age's value; developmental consequences in children are **not modelled** (flag `age_effect_unmodelled`).
- The person is assumed to start with the reference (average) alleles. Effects are applied to the average
  body, not to an individual's genome or ancestry.
- **Catalogue.** `data/genetics/catalogue.yaml`: MSTN loss of function (muscle), EPOR truncating gain of
  function (haemoglobin, haematocrit) and GHR loss of function (height). **All current effects are hypothetical:**
  the magnitudes are round-number guesses of a plausible order, from general knowledge that was not
  retrieved as a source in this project. Replace them with sourced values (tier *measured*, with source and
  evidence) as evidence is added. CYP2C9/ibuprofen (Phase 7) remains the only measured genotype effect.
- **Full gene list.** `tools/import_hgnc.py` turns an HGNC complete-set download into
  `data/genetics/gene_index.yaml` (valid symbols, no effects); once present, `--hypothesis` rejects unknown
  symbols. It was written without network access to genenames.org and tested only on a small sample file in
  the same layout, so run it on the real download and check the output.

Not modelled: delivery, editing efficiency (you supply the fraction), off-target or unintended edits,
immune responses, effects on traits that depend on many genes, interactions between edits, and any
consequence over time.

## Web UI

```
python -m humansim serve --open      # http://127.0.0.1:8765/
```

A single page with five views: **Body** (every parameter, editable; changed values highlighted; flags and a
"what changed and why" list), **Male vs female**, **Hormones** (cycle day slider with profile charts for
women; testosterone setting and the assumed feedback simulation for men), **Drugs** (IV concentration curves,
with CYP2C9 genotype for ibuprofen) and **Gene edits** (somatic or germline, catalogue, your own
hypotheses). The sidebar sets sex, age (0-20 y, 20 = adult) and whether edits ripple through the solver. A
yellow banner stays visible whenever the state contains hypothetical gene-edit effects.

It runs on your machine only: standard library, no extra dependencies, listens on 127.0.0.1, refuses
requests whose Host header is not loopback, stores nothing. Edits to values are in the parameter's unit
(shown in the table). The page works offline; the layout adapts to a phone-width screen and to light and dark mode.
Note that once any hypothetical edit is active, every value the solver moves is flagged
`downstream_of_hypothetical`, including values the edit may not really affect; the flag is conservative.

## Layout

```
humansim/
  data/values/       parameters and their values, by domain
  data/sources/      citations (sidecar; the model runs without it)
  data/validation/   published statistics that derived values must reproduce
  data/modifiers/    sex, hormone, genotype and drug effects (sex_female.yaml so far)
  registry.py        loading and validation
  derived.py         derived formulas
  context.py         Context, State and resolve()
  couplings.py       coupling records, deadband and the coupling functions
  solver.py          steady-state solve()
  dynamics.py        generic RK4 ODE engine
  hpg.py             male HPG feedback model (assumed structure)
  cycle.py           menstrual-cycle hormone profile
  pharmacology.py    one-compartment IV pharmacokinetics
  genetics.py        pharmacogene alleles, activity scores, frequencies
  editing.py         somatic and germline gene edits, measured or hypothetical
  server.py          local web server and JSON API
  web/index.html     the single-page UI
  data/genetics/     gene records (genes.yaml)
  data/drugs/        drug records (seed.yaml)
  data/couplings/    the Phase 3 couplings (steady_state.yaml)
  growth.py          CDC growth-curve (LMS) lookup for the age axis
  data/growth/       CDC 2000 LMS tables and curves.yaml
  modifiers.py       the modifier mechanism
  validation.py      validation-file loader
  checks.py          consistency and dataset checks
  report.py          coverage and state reports
tools/build_ansur_data.py   regenerates the ANSUR II data from the raw CSVs
tools/build_sex_rules.py    regenerates the sex rules from the repository data
tools/build_icrp_data.py    the ICRP 89 transcription, with its integrity checks
tests/
docs/SPEC.md
```

Regenerate the ANSUR II data (the CSVs are treated as untrusted input, so use `-I`):

```
python -I tools/build_ansur_data.py ANSUR_II_MALE_Public.csv ANSUR_II_FEMALE_Public.csv
```

## Licenses

Code: Apache-2.0 (`LICENSE`). Data under `humansim/data/`: CC BY 4.0 (see `DATA-LICENSE.md`).
Check each source's redistribution terms as you add data. This is not legal advice.
