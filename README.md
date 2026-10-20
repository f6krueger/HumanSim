# HumanSim

A Python library that models the human body as a network of interconnected anatomical and
physiological parameters, for population averages (not individuals). You set a context, edit any
value, and see what changes. See `docs/SPEC.md` for the full design.

**Status: Phase 1.** The core, derived layer, consistency checks, modifier machinery and the adult
male baseline exist. Sex rules, the solver, development, hormones, pharmacology and genetics are
specified but not built yet. The data is thin: run `python -m humansim coverage` to see exactly how
thin.

## Install and run

```
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest

python -m humansim show                                  # resolved reference male
python -m humansim show --set body_mass=90 --changed     # what changes if mass is 90 kg
python -m humansim show --set body_mass="198 lb" --changed
python -m humansim coverage --missing                    # what still needs a value
python -m humansim check                                 # consistency checks
```

```python
from humansim import Registry, Context, resolve

reg = Registry.load()
state = resolve(reg, Context(edits={"body_mass": 90.0}))
state.get("body_mass_index")      # derived, recomputed
state.flags["body_mass"]          # {'edited'}
```

## How values are marked

Every stored parameter has a status:

| Status | Meaning |
|---|---|
| `referenced` | Value read from a source that is named in `data/sources/sources.yaml` |
| `unsourced` | A value is present but I could not check it against a source |
| `placeholder` | A stand-in |
| `missing` | The parameter is defined but has no value yet |

`derived` values are computed from anchors and are never stored. They inherit `unsourced` and
`edited` flags from their inputs.

## What is real and what is not

- **Referenced:** height 176 cm, mass 73 kg, cardiac output 6.5 L/min, blood volume 5.3 L,
  haematocrit 0.45, hepatic flow fraction 0.255. These were read from pages that quote the ICRP 89
  reference adult male, **not from ICRP 89 itself**, which I could not open.
- **Unsourced:** organ and tissue masses (recalled as ICRP 89 values), plus typical figures such as
  core temperature, haemoglobin and SaO2. Resting metabolic rate is Mifflin-St Jeor at an assumed
  age of 30. All are marked `unsourced`. Verify them before relying on them.
- **Missing:** all body dimensions other than height and mass (arm span, girths, segment lengths and
  so on). They need a real survey source such as ANSUR II or NHANES.
- **Left out on purpose:** silk, rete, structural lipid and the cold and working energy figures in
  the original document, which describe a non-human body plan rather than a human reference.

## Layout

```
humansim/
  data/values/      parameters and their values, by domain
  data/sources/     citations (sidecar; the model runs without it)
  data/modifiers/   sex, hormone, genotype and drug effects (none shipped yet)
  registry.py       loading and validation
  derived.py        derived formulas
  context.py        Context, State and resolve()
  modifiers.py      the modifier mechanism
  checks.py         consistency and dataset checks
  report.py         coverage and state reports
tests/
docs/SPEC.md
```

## Licenses

Code: Apache-2.0 (`LICENSE`). Data under `humansim/data/`: CC BY 4.0 (see `DATA-LICENSE.md`).
Check each source's redistribution terms as you add data. This is not legal advice.
