"""Tight-binding models for group 4 transition metal dichalcogenides (tmd), 11 band."""
import numpy as np
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices


class TmdNN123MeoXeo(AbstractLattice):
    r"""Monolayer of a group 4 tmd using the second nearest-neighbor 11-band model

    Parameters
    ----------
    name : str
        Name of the tmd to model. The available options are: MoS2, WS2, MoSe2,
        WSe2. The relevant tight-binding parameters for these
        materials are given by https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    override_params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must
        be in the format `"name": [eps1,    eps3,    eps4,    eps6,    eps7,    eps9,   eps10,  t1_1_1,  t1_2_2,  t1_3_3,
                                   t1_4_4,  t1_5_5,  t1_6_6,  t1_7_7,  t1_8_8,  t1_9_9,t1_10_10,t1_11_11,  t1_3_5,  t1_6_8,
                                   t1_9_11,  t1_1_2,  t1_3_4,  t1_4_5,  t1_6_7,  t1_7_8, t1_9_10,t1_10_11,  t5_4_1,  t5_3_2,
                                   t5_5_2,  t5_9_6, t5_11_6, t5_10_7,  t5_9_8, t5_11_8,  t6_9_6, t6_11_6,  t6_9_8, t6_11_8,
                                   a,       c,     dXX,     dXM,    lambM, lambX]`.
    """
    def __init__(self, **kwargs):
        lattice_orbital_dict = {"l": {"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
                                "orbs": {"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"],
                                         "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
                                "group": {"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}}
        super().__init__(orbital=lattice_orbital_dict, n_v=6, n_b=11)
        self.lattice_name = "11 bands 3NN model"
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_0_x = self.block_diag(t_m.e_xe, t_m.e_xo)
        h_1_m = self.block_diag(t_m.t_1_me, t_m.t_1_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_2_x = self.block_diag(t_m.t_2_xe, t_m.t_2_xo)
        h_3_m = self.block_diag(t_m.t_3_me, t_m.t_1_mo * 0)
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "h_3_m", "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, h_3_m, self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))
