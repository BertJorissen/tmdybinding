"""Test the generation of the lattices."""
import pytest
import tmdybinding as tmdy
import pybinding as pb
import numpy as np

lattices = {
    "hexagonal-liu2-mos2": tmdy.TmdNN2Me(params=tmdy.liu2["MoS2"]).lattice(),
    "hexagonal-liu2-mose2": tmdy.TmdNN2Me(params=tmdy.liu2["MoSe2"]).lattice(),
    "hexagonal-liu2-mote2": tmdy.TmdNN2Me(params=tmdy.liu2["MoTe2"]).lattice(),
    "hexagonal-liu2-ws2": tmdy.TmdNN2Me(params=tmdy.liu2["WS2"]).lattice(),
    "hexagonal-liu2-wse2": tmdy.TmdNN2Me(params=tmdy.liu2["WSe2"]).lattice(),
    "hexagonal-liu2-wte2": tmdy.TmdNN2Me(params=tmdy.liu2["WTe2"]).lattice(),
    "hexagonal-liu6-mos2": tmdy.TmdNN2Me(params=tmdy.liu6["MoS2"]).lattice(),
    "hexagonal-liu6-mose2": tmdy.TmdNN2Me(params=tmdy.liu6["MoSe2"]).lattice(),
    "hexagonal-liu6-mote2": tmdy.TmdNN2Me(params=tmdy.liu6["MoTe2"]).lattice(),
    "hexagonal-liu6-ws2": tmdy.TmdNN2Me(params=tmdy.liu6["WS2"]).lattice(),
    "hexagonal-liu6-wse2": tmdy.TmdNN2Me(params=tmdy.liu6["WSe2"]).lattice(),
    "hexagonal-liu6-wte2": tmdy.TmdNN2Me(params=tmdy.liu6["WTe2"]).lattice(),
    "hexagonal-wu": tmdy.TmdNN256Meo(params=tmdy.wu["MoS2"]).lattice(),
    "hexagonal-fang-mos2": tmdy.TmdNN123MeoXeo(params=tmdy.fang["MoS2"]).lattice(),
    "hexagonal-fang-mose2": tmdy.TmdNN123MeoXeo(params=tmdy.fang["MoSe2"]).lattice(),
    "hexagonal-fang-ws2": tmdy.TmdNN123MeoXeo(params=tmdy.fang["WS2"]).lattice(),
    "hexagonal-fang-wse2": tmdy.TmdNN123MeoXeo(params=tmdy.fang["WSe2"]).lattice(),
    "hexagonal-jorissen": tmdy.TmdNN123MeoXeo(params=tmdy.jorissen["MoS2"]).lattice(),
    "hexagonal-all": tmdy.TmdNN2Me(params=tmdy.all["MoS2"]).lattice(),
    "hexagonal-rostami": tmdy.TmdNN12MeXe(params=tmdy.rostami["MoS2"]).lattice(),
    "hexagonal-dias-mos2": tmdy.TmdNN125MeoXeo(params=tmdy.dias["MoS2"]).lattice(),
    "hexagonal-dias-mose2": tmdy.TmdNN125MeoXeo(params=tmdy.dias["MoSe2"]).lattice(),
    "hexagonal-dias-mote2": tmdy.TmdNN125MeoXeo(params=tmdy.dias["MoTe2"]).lattice(),
    "hexagonal-dias-ws2": tmdy.TmdNN125MeoXeo(params=tmdy.dias["WS2"]).lattice(),
    "hexagonal-dias-wse2": tmdy.TmdNN125MeoXeo(params=tmdy.dias["WSe2"]).lattice(),
    "hexagonal-cappelluti": tmdy.TmdNN12MeoXeo(params=tmdy.cappelluti["MoS2"]).lattice(),
    "hexagonal-roldan": tmdy.TmdNN12MeoXeo(params=tmdy.roldan["MoS2"]).lattice()
}


def baseline_energy_calculator(lat4=False, soc=False, eo=False, so=False, sp=False):
    """Return the energies."""
    solver = pb.solver.lapack(pb.Model(
        tmdy.TmdNN123456MeoXeo(lat4=lat4, soc=soc, soc_eo_flip=eo, single_orbital=so, soc_polarized=sp).lattice(),
        pb.translational_symmetry(),
        pb.force_double_precision()
    ))
    solver.set_wave_vector([0.1, 0.2, 0.0])
    return np.sort(solver.eigenvalues).tolist()


energies = {
    "all": baseline_energy_calculator(),
    "all-lat4": baseline_energy_calculator(lat4=True),
    "all-soc": baseline_energy_calculator(soc=True),
    "all-soc-lat4": baseline_energy_calculator(soc=True, lat4=True),
    "all-soc-eo": baseline_energy_calculator(eo=True),
    "all-soc-eo-lat4": baseline_energy_calculator(eo=True, lat4=True),
    "all-sp": baseline_energy_calculator(sp=True),
    "all-lat4-sp": baseline_energy_calculator(lat4=True, sp=True),
    "all-soc-sp": baseline_energy_calculator(soc=True, sp=True),
    "all-soc-lat4-sp": baseline_energy_calculator(soc=True, lat4=True, sp=True),
    "all-soc-eo-sp": baseline_energy_calculator(soc=True, eo=True, sp=True),
    "all-soc-eo-lat4-sp": baseline_energy_calculator(soc=True, eo=True, lat4=True, sp=True)
}


@pytest.fixture(scope="module", ids=list(lattices.keys()), params=list(lattices.values()))
def lattice(request):
    return request.param


@pytest.fixture(scope="module", ids=list(energies.keys()), params=list(energies.values()))
def energy(request):
    return request.param


def test_expected(lattice, baseline):
    """Test the most basic lattice from Liu2."""
    expected = baseline(lattice)
    assert pytest.fuzzy_equal(lattice, expected)


def test_energies(energy, baseline):
    """Test the most basic lattice from Liu2."""
    expected = baseline(energy)
    assert pytest.fuzzy_equal(energy, expected)
