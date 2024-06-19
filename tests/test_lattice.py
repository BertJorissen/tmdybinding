"""Test the generation of the lattices."""
import pytest
import tmdybinding as tmdy

lattices = {
    "hexagonal-liu2": tmdy.TmdNN2Me(params=tmdy.liu2["MoSe2"]).lattice(),
    "hexagonal-all": tmdy.TmdNN2Me(params=tmdy.all["MoS2"]).lattice()
}


@pytest.fixture(scope="module", ids=list(lattices.keys()), params=list(lattices.values()))
def lattice(request):
    return request.param


def test_expected(lattice, baseline):
    """Test the most basic lattice from Liu2."""
    expected = baseline(lattice)
    assert pytest.fuzzy_equal(lattice, expected)
