# Couplings

`steady_state.yaml` holds the Phase 3 edges. See `humansim/couplings.py` for the model and
`README.md` for the summary. A record has `id`, `target` (or `targets`, or `select`), `fn`, `args`,
`kind` (mechanistic | empirical | assumed), `role` (accounting | compensation), `deadband`
(fraction or `ref_range`; compensation only), `validity` (range of the driving ratio), `status`,
`source_ref`, `note`. Each target may have at most one coupling. Values solved from a coupling
whose status is not `referenced` carry the flag `coupling_<status>`.
