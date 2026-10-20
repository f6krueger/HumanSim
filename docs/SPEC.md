# Human Reference Simulation: Specification (v0.3 draft)

Working name: `humansim` (rename any time). Code license: Apache-2.0. Data license: CC BY 4.0 (both easy to change; see section 14).

## 1. Purpose

A Python library that models the human body as a network of interconnected anatomical and physiological parameters, for population averages rather than individuals. A user sets a context (sex, age, hormone state, genotype, drugs), edits any value, and sees the consequences ripple through the system. Intended for researchers, clinicians, teachers, and anyone who wants a reference human to experiment on. A web UI comes after the library works.

## 2. Scope

All of the following are in scope for the design and the architecture from day one:

- **Sex:** adult reference male; female derived by dimorphism rules.
- **Development:** age axis from newborn to adult (and later, ageing).
- **Hormones:** concentrations and feedback axes that modulate physiology.
- **Pharmacology:** drug exposure and effect, built on the same organ masses and blood flows.
- **Genetics:** variants that shift parameters or drug handling.

**Delivery is staged so something useful exists first** (section 12). Each module ships with a small seed dataset, then grows as data is added. A module with thin data is still present and flagged, never omitted.

## 3. Design principles

1. **Real reference data is the baseline.** The supplied document supplies the parameter list and structure; its values are placeholders.
2. **One mechanism for all variation.** Sex, age, hormones, genotype, and drugs are all *modifiers* on a baseline (section 5). Adding a new kind of variation does not change the core.
3. **Derived values are never stored.** BMI, BSA, density, and indices are computed from anchors.
4. **Within normal variation, nothing changes.** Downstream adjustment starts only when a regulated parameter leaves its reference range.
5. **Say when extrapolating.** Output carries uncertainty and flags results outside the validity range of the data.
6. **Values and citations are separate.** The tool runs on values alone; sources are a sidecar dataset (section 8).

## 4. Parameter registry

Each parameter is one record in a values file:

| Field | Meaning |
|---|---|
| `id` | Stable machine name, e.g. `liver_mass` |
| `name`, `unit` | Display name; units handled with `pint` |
| `kind` | `anchor`, `derived`, `regulated`, or `dynamic` (time-varying, e.g. hormone level, drug concentration) |
| `value` | Mean, indexed by (sex, age) where available |
| `dispersion` | SD or percentile bounds, if known |
| `tiers` | Reference range, compensation range, viability limit |
| `source_ref` | Optional pointer into the sources dataset |
| `status` | `referenced`, `unsourced`, or `placeholder` (default `unsourced`) |

The `status` flag is the only citation information the working tool needs. Everything else lives in the sidecar.

## 5. Context and modifiers

The simulation takes a **context**:

```
Context(sex, age, hormone_state, genotype, drugs, edits)
```

Resolution order:
1. Look up the baseline for (sex, age) from the values.
2. Apply **modifiers** in a fixed order: sex rules, hormone effects, genotype effects, drug effects, then user edits.
3. Recompute derived values.
4. Solve the coupling graph.

A modifier is a small record: which parameters it touches, the form of the effect (ratio, offset, or function), its validity range, and its status. Examples:
- **Sex rule:** pelvic breadth ratio, blood volume ratio.
- **Hormone effect:** testosterone level shifts muscle mass and haemoglobin.
- **Genotype effect:** CYP2D6 poor-metaboliser status lowers clearance of its substrates.
- **Drug effect:** a drug changes heart rate or blood pressure through its dose-response.

## 6. Coupling graph and solver

- Nodes are parameters; edges are equations or rules, each labeled **empirical** (fitted, with validity range) or **mechanistic** (physical or physiological law).
- Engine: `networkx` for structure, `sympy` or plain callables for equations, `scipy` for constraint solving and ODEs.

**Two solve modes**, because some modules are time-dependent:
- **Steady-state:** for anatomy and physiology (organ sizes, flows, metabolic chain).
- **Dynamic:** ODE integration for hormones (feedback axes such as the hypothalamic-pituitary-gonadal axis), growth over age, and drug concentration over time.

**Three-tier ranges** for regulated parameters: reference range (deadband, nothing downstream changes), compensation zone (solver adjusts coupled parameters), viability limit (compensation fails, state flagged). Missing tier data returns "range undefined" instead of a guess.

**First chains:** (a) metabolic rate → oxygen demand → cardiac output → heart mass and blood volume → heat and water; (b) liver mass → hepatic blood flow → abdominal volume → lung compression → metabolism. The same organ masses and flows feed the pharmacology module, so (b) doubles as PBPK groundwork.

## 7. Sex derivation, development, hormones, pharmacology, genetics

- **Sex:** female = male baseline plus per-parameter dimorphism rules. Validated against the published reference female; mismatches trigger a rule fix.
- **Development:** values indexed by age; anchor ages from ICRP 89 (newborn, 1, 5, 10, 15, adult) plus WHO and CDC growth curves for continuous measures. Interpolated values are flagged. Pediatric and adult references are never silently mixed.
- **Hormones:** each hormone is a `dynamic` parameter with a reference range by sex, age, and (for menstrual cycle) phase. Effects on other parameters are modifiers. Feedback loops are ODE models, added axis by axis.
- **Pharmacology:** compartmental PBPK using the registry's organ masses, blood flows, and composition; drug records hold clearance, binding, and dose-response. Genetic modifiers adjust clearance.
- **Genetics:** a variant table (gene, allele, population frequency, effect on parameters or drug handling). Effects are modifiers; population frequencies let the tool produce population-weighted averages.

## 8. Citations as a separate dataset

- `data/values/`: the working data (parameters, modifiers). Minimal fields, loads fast, runs without citations.
- `data/sources/`: a separate dataset keyed by `source_ref`, holding full citation, edition, population, year, sample size, evidence grade (A to D), and notes.
- A **coverage report** command lists which values are `referenced`, `unsourced`, or `placeholder`, so you can see where citation work is needed.
- Sourcing can proceed in parallel with development: the tool is useful before it is fully cited, and citation work never blocks a code change.
- Corrections are made in the sources dataset and show up in output automatically.

## 9. Derived layer and consistency checks

Formulas for BMI, body surface area, density, segment fractions, limb indices, and so on. The checker verifies, for example, that tissue masses sum to body mass (reporting the unassigned remainder), that volume and density reconcile with mass, and that segments are consistent with stature.

## 10. Solver output

Per run: new values, change from baseline, uncertainty, and flags (outside data range, unsourced, viability exceeded, inconsistency, interpolated). The result also lists which modifiers and edges contributed to each changed value, so the ripple can be traced.

## 11. Repository layout (proposed)

```
humansim/
  data/
    values/        # parameters, modifiers (by domain)
    sources/       # citations and evidence grades
  core/            # registry, units, context, modifiers
  derived.py
  graph.py
  solver.py        # steady-state and dynamic
  modules/
    sex.py
    development.py
    hormones.py
    pharmacology.py
    genetics.py
  validate/
  report.py        # coverage and consistency reports
tests/
docs/
```

## 12. Build phases

Each phase ends with something usable.

1. **Core:** registry, units, context and modifier system, derived layer, real adult male baseline.
2. **Sex module:** both sexes stored; sex rules fill gaps in either direction; checked against an independent dataset. (Built.)
3. **Steady-state solver:** metabolic chain and liver test; deadband on compensation couplings, accounting couplings immediate, validity ranges. Built as reference-calibrated ratio couplings (see README). Three tiers in the original sense (reference / compensation / viability) are only partly present: viability limits are still unpopulated, so "range undefined" applies.
4. **Development module:** age axis with ICRP anchor ages and growth curves.
5. **Dynamic solver and hormones:** first axis (sex hormones), then others.
6. **Pharmacology:** PBPK on the existing anatomy, with a small seed drug set.
7. **Genetics:** variant table and its effects on parameters and drug clearance.
8. **Web UI.**
9. **Ongoing:** citation work, new data, new axes.

## 13. Risks

- **Scope is large.** The architecture handles it, but content volume does not shrink: every module needs real data. Staging and the coverage report are the control.
- **Dynamic models need care.** Hormone feedback and PBPK are well studied but parameter-hungry; early versions will be coarse and flagged as such.
- **Tier 2 and 3 ranges and many couplings are sparsely documented.** Expect "range undefined" often.
- **Extreme extrapolations** (a 9 kg liver) rest on mechanistic assumptions, not data.
- **Genetics adds a population layer.** Allele frequencies differ across populations, so the baseline must state which population it represents.
- **Clinical-adjacent use.** Add a clear notice that outputs are reference estimates, not medical advice, before any clinician-facing release.

## 14. Licensing notes

- Code: Apache-2.0 (permissive, with a patent grant). MIT is the simpler alternative.
- Data: CC BY 4.0 (attribution required) or CC0 (no restrictions). Individual facts are generally not copyrightable, but source tables and compilations may carry terms, so each source's redistribution terms should be checked as it is added. This is not legal advice.

## 15. Decisions for review

1. Apache-2.0 plus CC BY 4.0, or something different?
2. Is the staged order in section 12 right, or should a module move up (for example pharmacology before development)?
3. Anything to add to the scope list in section 2?
