"""Tight-binding models for group 1 transition metal dichalcogenides (tmd), 5 band."""
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices


class TmdNN256Me(AbstractLattice):
    r"""Monolayer of a group 6 tmd using the nearest-neighbor 3-band model

    Parameters
    ----------
    name : str
        Name of the tmd to model. The available options are: MoS2, WS2, MoSe2,
        WSe2, MoTe2, WTe2. The relevant tight-binding parameters for these 
        materials are given by https://doi.org/10.1103/PhysRevB.88.085433.
    override_params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must 
        be in the format `"name": [a, eps1, eps2, t0, t1, t2, t11, t12, t22]`.
    tnn: Boolean
        Take the Third Nearest Neighbor into account
    soc: Boolean
        Also calculate the Spin Orbit Coupling
    transform: (2x2) numpy array
        transformation matrix for the basis vectors
    """

    def __init__(self, **kwargs):

        lattice_orbital_dict = {"l": {"M": [0, 2, -2]},
                                "orbs": {"M": ["dz2", "dx2y2", "dxy"]},
                                "group": {"M": [0, 1, 1]}}
        super().__init__(orbital=lattice_orbital_dict, n_v=0, n_b=3)
        self.lattice_name = "3 bands 6NN model"
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = t_m.e_me
        h_2_m = t_m.t_2_me
        h_5_m = t_m.t_5_me
        h_6_m = t_m.t_6_me
        keys = ["h_0_m", "h_2_m", "h_5_m", "h_6_m", "a", "lamb_m"]
        values = [h_0_m, h_2_m, h_5_m, h_6_m, self.params["a"], self.params["lamb_m"]]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))
