"""Tight-binding models for transition metal dichalcogenides (tmd), 6 band."""
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices


class TmdNN12MeXe(AbstractLattice):
    r"""Monolayer of a group 4 tmd using the second nearest-neighbor 6-band model

    Parameters
    ----------
    name : str
        Name of the tmd to model. The available options are: MoS2, WS2, MoSe2,
        WSe2. The relevant tight-binding parameters for these 
        materials are given by https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    override_params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must 
        be in the format `
        "name": [
                  a,   eps_6,   eps_7,   eps_9, eps_10,  t1_6_6,   1_7_7,
             t1_8_8,  t1_9_9,t1_10_10,t1_11_11, t1_6_8, t1_9_11,   1_6_7,
             t1_7_8, t1_9_10,t1_10_11,  t5_9_6,t5_11_6, t5_10_7,  t5_9_8,
            t5_11_8,       c,    d_XX,    d_XM, lamb_m,  lamb_x
        ]`.
    """
    def __init__(self, **kwargs):

        lattice_orbital_dict = {"l": {"M": [0, 2, -2], "X": [1, -1, 0]},
                                "orbs": {"M": ["dz2", "dx2y2", "dxy"],
                                         "X": ["pxe", "pye", "pze"]},
                                "group": {"M": [0, 2, 2], "X": [0, 0, 1]}}
        super().__init__(orbital=lattice_orbital_dict, n_v=3, n_b=6)
        self.lattice_name = "6 bands model"
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = t_m.e_me
        h_0_x = t_m.e_xe
        h_1_m = t_m.t_1_me
        h_2_m = t_m.t_2_me
        h_2_x = t_m.t_2_xe
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))
