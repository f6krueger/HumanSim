"""Genetics (Phase 7): genotype to activity score, phenotype, population frequency and drug exposure.

Scope, plainly: pharmacogenetics for one gene (CYP2C9) feeding drug clearance, built on the CPIC
activity-score scheme. A person's genotype is `Context.genotype = {"CYP2C9": "*1/*3"}`. The model
does not change any anatomical or physiological parameter from genotype: no such effect with numbers
is in the registry's sources yet, and none is invented.

Reference convention: the reference person (the label values) is treated as a normal metabolizer
(*1/*1, activity score 2). A real population mean sits slightly above that, which is not corrected.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Mapping

import yaml

from .registry import DATA_DIR, RegistryError


@dataclass(frozen=True)
class Gene:
    id: str
    name: str
    alleles: Mapping[str, float]  # allele -> activity value
    functions: Mapping[str, str]
    phenotypes: Mapping[float, str]  # activity score -> phenotype
    populations: Mapping[str, Mapping]
    activity_source_ref: str

    def parse(self, diplotype: str) -> tuple[str, str]:
        parts = [a.strip() for a in str(diplotype).split("/")]
        if len(parts) != 2 or any(a not in self.alleles for a in parts):
            raise ValueError(f"{self.id}: diplotype must be two of {sorted(self.alleles)} like '*1/*3', got {diplotype!r}")
        return tuple(sorted(parts))  # type: ignore[return-value]

    def activity_score(self, diplotype: str) -> float:
        a, b = self.parse(diplotype)
        return self.alleles[a] + self.alleles[b]

    def phenotype(self, diplotype: str) -> str:
        score = self.activity_score(diplotype)
        return self.phenotypes[score]

    def allele_frequencies(self, population: str) -> dict[str, float]:
        if population not in self.populations:
            raise KeyError(f"{self.id}: no frequency data for population {population!r}; have {sorted(self.populations)}")
        given = dict(self.populations[population]["frequencies"])
        ref_allele = next(a for a, v in self.alleles.items() if v == max(self.alleles.values()))
        given[ref_allele] = given.get(ref_allele, 0.0) + 1.0 - sum(given.values())
        return given

    def diplotype_frequencies(self, population: str) -> dict[str, float]:
        """Hardy-Weinberg diplotype frequencies from the allele frequencies (assumes random mating)."""
        f = self.allele_frequencies(population)
        out = {}
        for a, b in combinations_with_replacement(sorted(f), 2):
            out[f"{a}/{b}"] = f[a] * f[b] * (1 if a == b else 2)
        return out

    def phenotype_frequencies(self, population: str) -> dict[str, float]:
        out: dict[str, float] = {}
        for dip, p in self.diplotype_frequencies(population).items():
            ph = self.phenotype(dip)
            out[ph] = out.get(ph, 0.0) + p
        return out


def load_genes(data_dir: str | Path | None = None) -> dict[str, Gene]:
    root = Path(data_dir) if data_dir else DATA_DIR
    out: dict[str, Gene] = {}
    for path in sorted((root / "genetics").glob("*.yaml")):
        for raw in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("genes") or []:
            alleles = {a: float(v["activity"]) for a, v in raw["alleles"].items()}
            g = Gene(id=raw["id"], name=raw["name"], alleles=alleles,
                     functions={a: v["function"] for a, v in raw["alleles"].items()},
                     phenotypes={float(k): v for k, v in raw["phenotypes"].items()},
                     populations=raw.get("populations") or {}, activity_source_ref=raw["activity_source_ref"])
            for dip in combinations_with_replacement(sorted(alleles), 2):
                if sum(alleles[a] for a in dip) not in g.phenotypes:
                    raise RegistryError(f"{g.id}: no phenotype for activity score of {'/'.join(dip)}")
            if g.id in out:
                raise RegistryError(f"Duplicate gene {g.id!r}")
            out[g.id] = g
    return out


def validate_genotype(genes: Mapping[str, Gene], genotype: Mapping[str, str]) -> dict[str, tuple[float, str]]:
    """Check a genotype and return {gene: (activity score, phenotype)}."""
    out = {}
    for gid, dip in genotype.items():
        if gid not in genes:
            raise KeyError(f"Unknown gene {gid!r}; known: {', '.join(genes)}")
        out[gid] = (genes[gid].activity_score(dip), genes[gid].phenotype(dip))
    return out


def exposure_ratio(table: Mapping[str, float], score: float) -> tuple[float, bool]:
    """AUC ratio at an activity score from a {score: ratio} table. Between two stored scores it is
    linear and returned as interpolated; outside the stored span there is no value."""
    pts = sorted((float(k), float(v)) for k, v in table.items())
    for s, v in pts:
        if s == score:
            return v, False
    lo = [(s, v) for s, v in pts if s < score]
    hi = [(s, v) for s, v in pts if s > score]
    if not lo or not hi:
        raise ValueError(f"no exposure data at activity score {score:g} (data cover {pts[0][0]:g}-{pts[-1][0]:g})")
    (s0, v0), (s1, v1) = lo[-1], hi[0]
    return v0 + (v1 - v0) * (score - s0) / (s1 - s0), True
