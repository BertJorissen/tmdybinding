"""Tight-binding models for group 1 transition metal dichalcogenides (tmd), 5 band."""
from .parameters import ParametersList
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices

_default_5band_params = ParametersList(dict(zip(
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_5_0_m_o",    "u_5_2_m_o",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e",
        "u_6_0_m_o",    "u_6_1_m_o",    "u_6_2_m_o",
    ],
    [
        0.3190,         0.073,          None,           "MoS2",
        0.683,          1.707,          3.558,
        -0.146,         0.506,          0.114,          0.073,          0.162,          0.085,
        -0.189,         -0.024,         -0.117,
        0.060,          0.273,          -0.034,         0.167,          -0.077,
        -0.063,         0.025,
        -0.038,         0.001,          -0.046,          -0.150,         -0.176,         0.266,
        0.165,          0.140,          -0.122
    ]
)))


class Group1Tmd5Band(AbstractLattice):
    r"""Monolayer of a group 1 tmd using the third nearest-neighbor 5-band model

    Parameters
    ----------
    params : ParametersList
        Parameters for the model.
    override_params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must 
        be in the format `"name": [a, eps1, eps2, t0, t1, t2, t11, t12, t22]`.
    soc: Boolean
        Also calculate the Spin Orbit Coupling
    transform: (2x2) numpy array
        transformation matrix for the basis vectors
    """
    def __init__(self, **kwargs):

        lattice_orbital_dict = {"l": {"M": [0, 2, -2, 1, -1]},
                                "orbs": {"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"]}}
        super().__init__(orbital=lattice_orbital_dict, n_v=0, n_b=5)
        self.single_orbital = False
        self._berry_phase_factor = 1
        self.params = _default_5band_params
        self.lattice_name = "Liu/Wu 5 bands 6NN model"
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_5_m = self.block_diag(t_m.t_5_me, t_m.t_5_mo)
        h_6_m = self.block_diag(t_m.t_6_me, t_m.t_6_mo)
        keys = ["h_0_m", "h_2_m", "h_5_m", "h_6_m", "a", "lamb_m"]
        values = [h_0_m, h_2_m, h_5_m, h_6_m, self.params["a"], self.params["lamb_m"]]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))
