"""Tight-binding models for group 5 transition metal dichalcogenides (tmd), 11 band."""
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices


class TmdNN125MeoXeo(AbstractLattice):
    r"""Monolayer of a group 1 tmd using the second nearest-neighbor 11-band model

    Parameters
    ----------
    name : str
        Name of the tmd to model. The available options are: MoS2
        WThe relevant tight-binding parameters for these
        materials are given by https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.195402
    override_params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must be
        in the format `"name": [ eps0,   eps2,    epsp,     epsz,   V_pdd,   V_pdp,   V_ddd,   V_ddp,  V_ddde,  V_ppd,
                                V_ppp,      a,   lamXX,    lamXM,   lamMM,    lamM,    lamX]
    sz : float
        Z-component of the spin degree of freedom
    """
    def __init__(self, **kwargs):

        lattice_orbital_dict = {"l": {"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
                                "orbs": {"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"],
                                         "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
                                "group": {"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}}
        super().__init__(orbital=lattice_orbital_dict, n_v=6, n_b=11)
        self.lattice_name = "SK 11 bands 5NN model"
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_0_x = self.block_diag(t_m.e_xe, t_m.e_xo)
        h_1_m = self.block_diag(t_m.t_1_me, t_m.t_1_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_2_x = self.block_diag(t_m.t_2_xe, t_m.t_2_xo)
        h_5_m = self.block_diag(t_m.t_5_me, t_m.t_5_mo)
        h_5_x = self.block_diag(t_m.t_5_xe, t_m.t_5_xo)
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "h_5_m", "h_5_c", "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, h_5_m, h_5_x, self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))
