import pytest
import yaml

from humansim import Registry, RegistryError
from humansim.checks import dataset_checks
from humansim.units import UnitError, convert, parse_value, validate_unit


def test_registry_loads_and_ids_are_unique(registry):
    ids = registry.ids()
    assert len(ids) == len(set(ids))
    assert len(ids) > 100


def test_every_unit_parses(registry):
    for p in registry:
        validate_unit(p.unit)


def test_status_matches_values(registry):
    for p in registry:
        if p.kind == "derived":
            assert p.status == "derived" and not p.values
        else:
            assert (p.status == "missing") == (not p.values), p.id


def test_referenced_values_have_existing_sources(registry, sources):
    for r in dataset_checks(registry, sources):
        assert r.status != "fail", r.message


def test_user_parameter_list_is_present(registry):
    # A sample of ids from the original parameter document must exist.
    for pid in ("height", "body_mass", "heart_mass", "brain_mass", "blood_volume",
                "hemoglobin", "femoral_shaft_diameter", "molar_bite_force", "arm_span"):
        assert pid in registry


def test_reference_male_anchor_values(registry):
    # Height and mass come from ANSUR II (computed from the raw file); the cardiovascular
    # values are the ICRP 89 reference male, read via secondary papers.
    assert registry["height"].baseline("male", "adult") == pytest.approx(175.6, abs=0.3)
    assert registry["body_mass"].baseline("male", "adult") == pytest.approx(85.5, abs=0.5)
    assert registry["cardiac_output"].baseline("male", "adult") == 6.5
    assert registry["blood_volume"].baseline("male", "adult") == 5.3


def test_female_data_stored_where_available(registry):
    assert registry["height"].baseline("female", "adult") is not None
    assert registry["brain_metabolic_load"].baseline("female", "adult") is None


def test_loader_rejects_bad_files(tmp_path):
    (tmp_path / "values").mkdir()

    def write(params):
        (tmp_path / "values" / "x.yaml").write_text(yaml.safe_dump({"group": "g", "parameters": params}))

    write([{"id": "a", "name": "A", "unit": "kg", "kind": "anchor", "typo_key": 1}])
    with pytest.raises(RegistryError, match="unknown keys"):
        Registry.load(tmp_path)

    write([{"id": "a", "name": "A", "unit": "furlongs-per-banana", "kind": "anchor"}])
    with pytest.raises(RegistryError, match="Unknown unit"):
        Registry.load(tmp_path)

    write([{"id": "a", "name": "A", "unit": "kg", "kind": "anchor", "status": "missing",
            "values": {"male": {"adult": 1.0}}}])
    with pytest.raises(RegistryError, match="missing"):
        Registry.load(tmp_path)

    write([{"id": "a", "name": "A", "unit": "kg", "kind": "anchor"},
           {"id": "a", "name": "A2", "unit": "kg", "kind": "anchor"}])
    with pytest.raises(RegistryError, match="duplicate"):
        Registry.load(tmp_path)


def test_unit_helpers():
    assert convert(198, "lb", "kg") == pytest.approx(89.81, abs=0.01)
    assert parse_value("90", "kg") == 90.0
    assert parse_value("198 lb", "kg") == pytest.approx(89.81, abs=0.01)
    assert convert(37, "degC", "degF") == pytest.approx(98.6)
    with pytest.raises(UnitError):
        parse_value("5 seconds", "kg")
