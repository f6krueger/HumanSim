import pytest

from humansim import Context, Registry, load_modifiers, load_sources, resolve
from humansim.genetics import exposure_ratio, load_genes, validate_genotype
from humansim.pharmacology import Dose, genotype_clearance_factor, load_drugs, simulate
from humansim.registry import RegistryError


@pytest.fixture(scope="module")
def env():
    return Registry.load(), load_modifiers(), load_genes(), load_drugs()


def test_activity_scores_and_phenotypes_follow_cpic(env):
    g = env[2]["CYP2C9"]
    expect = {"*1/*1": (2, "normal"), "*1/*2": (1.5, "intermediate"), "*1/*3": (1, "intermediate"),
              "*2/*2": (1, "intermediate"), "*2/*3": (0.5, "poor"), "*3/*3": (0, "poor")}
    for dip, (score, ph) in expect.items():
        assert g.activity_score(dip) == score and g.phenotype(dip) == ph
    assert g.activity_score("*3/*1") == 2 - 1  # order does not matter


def test_bad_genotypes_are_refused(env):
    genes = env[2]
    with pytest.raises(KeyError):
        validate_genotype(genes, {"NOPE": "*1/*1"})
    for bad in ("*1", "*1/*9", "*1/*2/*3", ""):
        with pytest.raises(ValueError):
            genes["CYP2C9"].parse(bad)


def test_population_frequencies_are_hardy_weinberg(env):
    g = env[2]["CYP2C9"]
    f = g.allele_frequencies("european")
    assert f == pytest.approx({"*1": 0.82, "*2": 0.11, "*3": 0.07})
    d = g.diplotype_frequencies("european")
    assert sum(d.values()) == pytest.approx(1.0)
    assert d["*1/*1"] == pytest.approx(0.82 ** 2) and d["*1/*3"] == pytest.approx(2 * 0.82 * 0.07)
    assert sum(g.phenotype_frequencies("european").values()) == pytest.approx(1.0)
    with pytest.raises(KeyError):
        g.allele_frequencies("martian")


def test_exposure_ratio_interpolates_and_refuses_extrapolation():
    table = {"2": 1.0, "1": 1.43}
    assert exposure_ratio(table, 1) == (1.43, False)
    r, interp = exposure_ratio(table, 1.5)
    assert r == pytest.approx(1.215) and interp
    with pytest.raises(ValueError):
        exposure_ratio(table, 0)


def test_genotype_changes_ibuprofen_clearance_by_the_cpic_ratio(env):
    reg, mods, genes, drugs = env
    base = simulate(reg, resolve(reg, Context(), mods), drugs["ibuprofen"], [Dose(800, 0, 1)], 48, mods)
    var = simulate(reg, resolve(reg, Context(genotype={"CYP2C9": "*1/*3"}), mods), drugs["ibuprofen"],
                   [Dose(800, 0, 1)], 48, mods)
    assert var.clearance_l_h == pytest.approx(base.clearance_l_h / 1.43)
    assert var.auc_mg_h_l() == pytest.approx(base.auc_mg_h_l() * 1.43, rel=0.01)
    assert "genotype_effect" in var.flags and "genotype_effect" not in base.flags
    assert base.auc_mg_h_l() == pytest.approx(192.8, rel=0.02)  # reproduces the label AUC for 800 mg


def test_poor_metabolizer_has_no_ibuprofen_data_and_is_refused(env):
    reg, mods, genes, drugs = env
    with pytest.raises(ValueError):
        simulate(reg, resolve(reg, Context(genotype={"CYP2C9": "*3/*3"}), mods), drugs["ibuprofen"],
                 [Dose(800, 0, 1)], 24, mods)


def test_drug_without_gene_data_is_unchanged_but_noted(env):
    reg, mods, genes, drugs = env
    f, notes, flags = genotype_clearance_factor(drugs["midazolam"], {"CYP2C9": "*1/*3"}, genes)
    assert f == 1.0 and not flags and "no data" in notes[0]


def test_genetics_data_is_consistent(env):
    genes, drugs = env[2], env[3]
    src = load_sources()
    for g in genes.values():
        assert g.activity_source_ref in src
        assert all(p["source_ref"] in src for p in g.populations.values())
    for d in drugs.values():
        for e in d.genotype_effects:
            assert e["gene"] in genes and e["source_ref"] in src
    bad = """
genes:
  - {id: X, name: X, alleles: {"*1": {activity: 1, function: normal}, "*2": {activity: 0.25, function: decreased}},
     activity_source_ref: s, phenotypes: {"2": normal}}
"""
    import tempfile, pathlib
    from humansim.genetics import load_genes as lg
    with tempfile.TemporaryDirectory() as d:
        (pathlib.Path(d) / "genetics").mkdir()
        (pathlib.Path(d) / "genetics" / "x.yaml").write_text(bad)
        with pytest.raises(RegistryError):
            lg(d)
