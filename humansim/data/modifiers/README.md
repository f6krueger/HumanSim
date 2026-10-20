# Modifiers

Phase 1 ships no modifiers. Phase 2 (sex rules) adds the first ones here.

A modifier file looks like this:

```yaml
modifiers:
  - id: example_sex_ratio
    stage: sex            # sex | hormone | genotype | drug
    op: ratio             # ratio | offset | set
    effects: {body_mass: 0.82}
    when: {sex: female}
    status: unsourced     # referenced | unsourced | placeholder
    source_ref: null
    note: ""
```
