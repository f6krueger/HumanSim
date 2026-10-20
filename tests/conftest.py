import pytest

from humansim import Registry, load_sources


@pytest.fixture(scope="session")
def registry():
    return Registry.load()


@pytest.fixture(scope="session")
def sources():
    return load_sources()
