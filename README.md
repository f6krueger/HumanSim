# HumanSim

A Python library that models the human body as a network of interconnected anatomical and
physiological parameters, for population averages (not individuals). You set a context, edit any
value, and see what changes. See `docs/SPEC.md` for the full design.

**Status: Phase 1.** The core, derived layer, consistency checks, modifier machinery and the adult
male baseline exist. Sex rules, the solver, development, hormones, pharmacology and genetics are
specified but not built yet. Run `python -m humansim coverage` to see exactly how much data there is.

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

- **ANSUR II (computed from the raw files, 32 parameters):** height, body mass, and the lengths,
  breadths and girths that ANSUR II measures, with SD and 2.5-97.5 percentile range.
  `tools/build_ansur_data.py` regenerates `data/values/anthropometry.yaml` from the public CSVs.
  **Caveat: the sample is US Army personnel (n=4,082 men, mean age 30), not the general public.**
  Their mean mass is 85.5 kg.
- **ICRP 89, read via secondary papers (4 parameters):** cardiac output 6.5 L/min, blood volume
  5.3 L, haematocrit 0.45, hepatic flow fraction 0.255. The ICRP 89 original was not read. These
  describe a 73 kg reference man.
- **Unsourced (19):** organ and tissue masses (recalled as ICRP 89 values, also for a 73 kg man)
  and typical figures such as core temperature and haemoglobin. Verify before relying on them.
- **Missing (24):** no ANSUR II equivalent (femur and tibia length, knee girth, standing reach,
  head height, neck length, body density...) or no source yet (bone density, enamel, bite force...).
- **Known inconsistency:** the organ and cardiovascular data are for a 73 kg man while the
  dimensions are for an 85.5 kg man. `python -m humansim check` reports this as a warning. Reconciling
  them (scaling organ masses and flows to body mass) is the job of the Phase 3 solver.
- **Left out on purpose:** silk, rete, structural lipid and the cold and working energy figures in
  the original document, which describe a non-human body plan.

## Female data

`data/validation/ansur2_female.yaml` holds ANSUR II female statistics (n=1,986) and each ratio to
the male mean. They are **not** baseline values. In Phase 2 the female is derived from the male by
sex rules, and this file is what the derivation must reproduce.

## Layout

```
humansim/
  data/values/       parameters and their values, by domain
  data/sources/      citations (sidecar; the model runs without it)
  data/validation/   published statistics that derived values must reproduce
  data/modifiers/    sex, hormone, genotype and drug effects (none shipped yet)
  registry.py        loading and validation
  derived.py         derived formulas
  context.py         Context, State and resolve()
  modifiers.py       the modifier mechanism
  validation.py      validation-file loader
  checks.py          consistency and dataset checks
  report.py          coverage and state reports
tools/build_ansur_data.py   regenerates the ANSUR II data
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
