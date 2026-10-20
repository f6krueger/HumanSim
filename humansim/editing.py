"""Gene edits (somatic in vivo and germline), measured or hypothetical.

An edit changes parameters through the same modifier mechanism as everything else (stage "genotype"),
so it flows through the solver and into pharmacology. Every effect carries a tier:

  measured      a sourced effect size (source, evidence grade and the people it was measured in)
  hypothetical  NOT data. Either a catalogue entry whose magnitude is a stated guess, or an effect the
                user typed in as their own hypothesis. Flagged `hypothetical` on every value it moves and
                on everything the solver then moves because of it.

Nothing here predicts what an edit would do to a real person. The model describes an average adult
(or child) of the chosen sex, assumes the person starts with the reference (average) alleles, and
applies each effect as a ratio or offset on the average value.

Edit modes:
  somatic   in vivo editing of a fraction of cells. The effect is scaled by the edited fraction. That
            linear scaling is an ASSUMPTION (a locally edited tissue acting on a whole-body value in
            proportion to the share of cells edited); flagged `mosaic_assumed_linear`.
  germline  every cell carries the edit from conception (fraction 1). Heritable and, in many countries,
            prohibited or tightly restricted for human embryos. The model applies the same effect at
            any age; developmental effects in children are not modelled (flagged `age_effect_unmodelled`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import yaml

from .modifiers import Modifier
from .registry import DATA_DIR, Registry, RegistryError

TIERS = ("measured", "hypothetical")
MODES = ("somatic", "germline")
EDIT_OPS = ("ratio", "offset", "set")


@dataclass(frozen=True)
class Effect:
    param: str
    op: str
    value: float
    tier: str = "hypothetical"
    source_ref: str | None = None
    evidence: str | None = None
    ancestry: str | None = None
    rationale: str = ""

    def __post_init__(self) -> None:
        if self.tier not in TIERS:
            raise RegistryError(f"Effect on {self.param!r}: tier must be one of {TIERS}")
        if self.op not in EDIT_OPS:
            raise RegistryError(f"Effect on {self.param!r}: op must be one of {EDIT_OPS}")
        if self.tier == "measured" and not (self.source_ref and self.evidence):
            raise RegistryError(f"Effect on {self.param!r}: a measured effect needs source_ref and evidence")
        if self.tier == "hypothetical" and not self.rationale:
            raise RegistryError(f"Effect on {self.param!r}: a hypothetical effect needs a rationale (why this size)")
        if self.op == "ratio" and self.value <= 0:
            raise RegistryError(f"Effect on {self.param!r}: a ratio must be positive")


@dataclass(frozen=True)
class Locus:
    id: str
    gene: str
    name: str
    states: Mapping[str, tuple[str, tuple[Effect, ...]]]  # state -> (description, effects)


@dataclass(frozen=True)
class GeneEdit:
    """One edit. Give `locus` and `state` for a catalogue entry, or `hypothesis` (a tuple of Effect) for
    your own idea, or both (their effects are added together as separate modifiers)."""
    gene: str
    locus: str | None = None
    state: str | None = None
    hypothesis: tuple[Effect, ...] = ()
    mode: str = "somatic"
    fraction: float = 1.0
    tissue: str | None = None

    def __post_init__(self) -> None:
        if self.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}")
        if not 0 < self.fraction <= 1:
            raise ValueError("edited fraction must be above 0 and at most 1")
        if self.mode == "germline" and self.fraction != 1:
            raise ValueError("a germline edit is in every cell: fraction must be 1")
        if bool(self.locus) != bool(self.state):
            raise ValueError("give both locus and state, or neither")
        if not self.locus and not self.hypothesis:
            raise ValueError("an edit needs a catalogue locus/state or a hypothesis")


def load_catalogue(data_dir: str | Path | None = None) -> dict[str, Locus]:
    root = Path(data_dir) if data_dir else DATA_DIR
    out: dict[str, Locus] = {}
    for path in sorted((root / "genetics").glob("catalogue*.yaml")):
        for raw in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("loci") or []:
            states = {}
            for sname, s in (raw.get("states") or {}).items():
                effects = tuple(Effect(**e) for e in s.get("effects") or [])
                states[sname] = (s.get("description", ""), effects)
            loc = Locus(id=raw["id"], gene=raw["gene"], name=raw.get("name", raw["gene"]), states=states)
            if loc.id in out:
                raise RegistryError(f"Duplicate locus {loc.id!r}")
            out[loc.id] = loc
    return out


def load_gene_index(data_dir: str | Path | None = None) -> dict[str, dict]:
    """Optional index of known gene symbols (see tools/import_hgnc.py). Empty if not imported."""
    root = Path(data_dir) if data_dir else DATA_DIR
    path = root / "genetics" / "gene_index.yaml"
    if not path.exists():
        return {}
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return {g["symbol"]: g for g in doc.get("genes") or []}


def catalogue_checks(registry: Registry, catalogue: Mapping[str, Locus], sources: Mapping[str, dict]) -> list[str]:
    """Problems in the catalogue (empty list = fine)."""
    problems = []
    for loc in catalogue.values():
        for sname, (_, effects) in loc.states.items():
            for e in effects:
                where = f"{loc.id}/{sname}/{e.param}"
                if e.param not in registry:
                    problems.append(f"{where}: unknown parameter")
                elif registry[e.param].kind == "derived":
                    problems.append(f"{where}: derived parameter (edit an input)")
                if e.tier == "measured" and e.source_ref not in sources:
                    problems.append(f"{where}: unknown source {e.source_ref!r}")
    return problems


def _scaled(effect: Effect, fraction: float) -> float:
    if effect.op == "ratio":
        return 1.0 + fraction * (effect.value - 1.0)
    if effect.op == "offset":
        return fraction * effect.value
    if fraction != 1.0:
        raise ValueError(f"{effect.param}: a 'set' effect cannot be applied to a fraction of cells; use ratio or offset")
    return effect.value


def edit_modifiers(edits, catalogue: Mapping[str, Locus] | None = None, gene_index: Mapping | None = None) -> list[Modifier]:
    """Turn gene edits into genotype-stage modifiers."""
    catalogue = load_catalogue() if catalogue is None else catalogue
    gene_index = load_gene_index() if gene_index is None else gene_index
    mods: list[Modifier] = []
    for n, ed in enumerate(edits):
        if gene_index and ed.gene not in gene_index:
            raise KeyError(f"Unknown gene symbol {ed.gene!r} (not in the gene index)")
        effects: list[tuple[Effect, str]] = []
        if ed.locus:
            if ed.locus not in catalogue:
                raise KeyError(f"Unknown locus {ed.locus!r}; known: {', '.join(catalogue)}")
            loc = catalogue[ed.locus]
            if ed.state not in loc.states:
                raise KeyError(f"{ed.locus}: unknown state {ed.state!r}; known: {', '.join(loc.states)}")
            effects += [(e, f"catalogue {ed.locus}/{ed.state}") for e in loc.states[ed.state][1]]
        effects += [(e, "user hypothesis") for e in ed.hypothesis]
        for i, (e, origin) in enumerate(effects):
            where = f"{ed.mode}" + (f", {ed.fraction:.0%} of cells" if ed.mode == "somatic" else "")
            tag = ("MEASURED " + f"({e.source_ref}, evidence {e.evidence}" + (f", {e.ancestry}" if e.ancestry else "") + ")"
                   if e.tier == "measured" else "HYPOTHETICAL, not data")
            note = f"gene edit {ed.gene} [{origin}; {where}]: {tag}. {e.rationale}".strip()
            extra = []
            if ed.mode == "somatic" and ed.fraction < 1:
                extra.append("mosaic_assumed_linear")
            if ed.mode == "germline":
                extra.append("germline")
            mods.append(Modifier(id=f"edit:{ed.gene}:{n}.{i}", stage="genotype", op=e.op,
                                 effects={e.param: _scaled(e, ed.fraction)}, when={},
                                 status="referenced" if e.tier == "measured" else "hypothetical",
                                 source_ref=e.source_ref, note=note + (" Flags: " + ", ".join(extra) if extra else ""),
                                 rule_type=""))
    return mods
