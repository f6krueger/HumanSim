import pytest

from humansim import Context, load_modifiers, resolve

AGES = [0, 1, 5, 10, 15, "adult"]


def v(registry, pid, sex, age):
    return registry[pid].baseline(sex, age)


def test_values_checked_against_the_pdf_page(registry):
    # read off page 18 of the ICRP 89 PDF (Table 2.8) and checked against the page image
    assert v(registry, "adipose_tissue_mass", "female", "adult") == 22.5  # not the 18.0 first recalled
    assert v(registry, "skin_mass", "female", "adult") == 2.3  # not the 1.7 first recalled
    assert v(registry, "liver_mass", "male", "adult") == 1.8 and v(registry, "liver_mass", "female", "adult") == 1.4
    assert v(registry, "brain_mass", "male", 5) == 1.31 and v(registry, "brain_mass", "female", 5) == 1.18  # "1310/1180"
    assert v(registry, "lungs_mass", "male", "adult") == 0.5  # tissue only; with blood is 1.2
    assert v(registry, "lungs_mass_with_blood", "male", "adult") == 1.2
    assert v(registry, "heart_mass", "male", "adult") == 0.33


def test_table_2_8_columns_sum_to_the_reports_96_percent(registry):
    from humansim.checks import SEX_ORGANS, WHOLE_BODY_ROWS
    for sex in ("male", "female"):
        for age in AGES:
            total = sum(v(registry, i, sex, age) or 0 for i in WHOLE_BODY_ROWS)
            mass = v(registry, "reference_body_mass", sex, age)
            assert 0.94 < total / mass < 0.97, (sex, age, total / mass)


def test_skeleton_components_add_up(registry):
    parts = ["cortical_bone_mass", "trabecular_bone_mass", "active_marrow_mass", "inactive_marrow_mass",
             "cartilage_mass", "teeth_mass", "skeleton_miscellaneous_mass"]
    for sex in ("male", "female"):
        for age in AGES:
            assert sum(v(registry, p, sex, age) for p in parts) == pytest.approx(v(registry, "skeleton_mass", sex, age), rel=0.03)


def test_regional_flows_and_blood_contents_sum_to_100_percent(registry):
    flows = [p.id for p in registry if p.id.startswith("blood_flow_fraction_")]
    skeleton_parts = {"red_marrow", "trabecular_bone", "cortical_bone", "other_skeleton"}
    contents = [p.id for p in registry if p.id.startswith("blood_content_fraction_")]
    for sex in ("male", "female"):
        # liver counted as arterial flow only: the portal flow is the gut's blood passing through again
        flow = sum(v(registry, i, sex, "adult") for i in flows if i.rsplit("fraction_", 1)[1] not in skeleton_parts)
        flow += v(registry, "hepatic_arterial_flow_fraction", sex, "adult")
        assert flow == pytest.approx(1.0, abs=1e-3), sex
        assert sum(v(registry, i, sex, "adult") for i in contents) == pytest.approx(1.0, abs=5e-4), sex


def test_female_hepatic_flow_is_27_percent_not_the_male_value(registry):
    assert v(registry, "hepatic_flow_fraction", "male", "adult") == 0.255
    assert v(registry, "hepatic_flow_fraction", "female", "adult") == 0.27


def test_physiology_by_age(registry):
    assert [v(registry, "cardiac_output", "male", a) for a in AGES] == [0.6, 1.2, 3.4, 5.0, 6.1, 6.5]
    assert v(registry, "cardiac_output", "female", "adult") == 5.9 and v(registry, "cardiac_output", "female", 15) == 6.1
    assert v(registry, "resting_metabolic_rate", "male", "adult") == pytest.approx(68 * 24)
    assert v(registry, "resting_metabolic_rate", "female", "adult") == pytest.approx(52 * 24)
    assert v(registry, "blood_volume", "female", "adult") == 3.9 and v(registry, "blood_volume", "male", 5) == 1.4
    assert v(registry, "total_lung_capacity", "male", "adult") == 7.0 and v(registry, "total_lung_capacity", "male", 0.25) == 0.28
    assert v(registry, "water_requirement", "male", "adult") == 2.6
    assert v(registry, "cortical_bone_density", "male", 0) == 1.65


def test_all_icrp_values_are_referenced_to_icrp89(registry):
    icrp = [p for p in registry if p.source_ref == "icrp89"]
    assert len(icrp) > 130
    assert all(p.status == "referenced" for p in icrp)


def test_metabolic_chain_is_consistent_with_organ_data(registry):
    # the organ-sum metabolic rate (Wang 2010 rates) vs the ICRP basal rate, two independent sources
    from humansim.checks import rmr_organ_crosscheck, mifflin_crosscheck
    for sex in ("male", "female"):
        s = resolve(registry, Context(sex=sex), load_modifiers())
        assert abs(rmr_organ_crosscheck(s).details["difference"]) < 0.10
        assert mifflin_crosscheck(s).status == "pass"


def test_red_cell_volume_cross_check(registry):
    from humansim.checks import icrp_red_cell_crosscheck
    m = icrp_red_cell_crosscheck(resolve(registry, Context(), load_modifiers()))
    assert m.status == "pass"
    f = icrp_red_cell_crosscheck(resolve(registry, Context(sex="female"), load_modifiers()))
    # 3.9 L x 0.42 = 1.64 L against ICRP's 1.5 L: the two sources disagree by about 9%, and it is reported
    assert f.details["model"] / f.details["icrp"] - 1 == pytest.approx(0.093, abs=0.01)
