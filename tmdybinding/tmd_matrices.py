""" The matrices for TMD 1H-phases and symmetry arguments

Date: 10 August 2023
Author: Bert Jorissen
"""

import numpy as np
from .sg_parameters import ParametersList


class TmdMatrices:
    """Construct the TMD hopping matrices"""
    def __init__(self, params: ParametersList):
        self.params = params

    @staticmethod
    def _t_n_e(u_0, u_1, u_2, u_3, u_4):
        return np.array([[0, 0, u_0], [u_1, u_2, 0], [u_3, u_4, 0]])

    @staticmethod
    def _t_n_o(u_0, u_1, u_2):
        return np.array([[u_0, 0], [0, u_1], [0, u_2]])

    @staticmethod
    def _t_m_xr(u_0, u_1, u_2, u_3, u_4, u_5):
        return np.array([[u_0, u_1, u_2], [-u_1, u_3, u_4], [-u_2, u_4, u_5]])

    @staticmethod
    def _t_m_me(u_0, u_1, u_2, u_3, u_4, u_5):
        return np.array([[u_0, u_1, u_2], [u_1, u_3, u_4], [-u_2, -u_4, u_5]])

    @staticmethod
    def _t_m_mo(u_0, u_1, u_2):
        return np.array([[u_0, u_1], [-u_1, u_2]])

    @staticmethod
    def _t_5_xr(u_0, u_2, u_3, u_5, u_6):
        return np.array([[u_3, 0, 0], [0, u_0, u_2], [0, u_6, u_5]])

    @property
    def e_xe(self):
        return np.diag((self.params["eps_0_x_e"], self.params["eps_0_x_e"], self.params["eps_1_x_e"]))

    @property
    def e_xo(self):
        return np.diag((self.params["eps_0_x_o"], self.params["eps_0_x_o"], self.params["eps_1_x_o"]))

    @property
    def e_me(self):
        return np.diag((self.params["eps_0_m_e"], self.params["eps_1_m_e"], self.params["eps_1_m_e"]))

    @property
    def e_mo(self):
        return np.diag((self.params["eps_0_m_o"], self.params["eps_0_m_o"]))

    @property
    def t_1_me(self):
        return self._t_n_e(*[self.params[f"u_1_{i}_m_e"] for i in range(5)])

    @property
    def t_1_mo(self):
        return self._t_n_o(*[self.params[f"u_1_{i}_m_o"] for i in range(3)])

    @property
    def t_2_me(self):
        return self._t_m_me(*[self.params[f"u_2_{i}_m_e"] for i in range(6)])

    @property
    def t_2_mo(self):
        return self._t_m_mo(*[self.params[f"u_2_{i}_m_o"] for i in range(3)])

    @property
    def t_2_xe(self):
        return self._t_m_xr(*[self.params[f"u_2_{i}_x_e"] for i in range(6)])

    @property
    def t_2_xo(self):
        return self._t_m_xr(*[self.params[f"u_2_{i}_x_o"] for i in range(6)])

    @property
    def t_3_me(self):
        return self._t_n_e(*[self.params[f"u_3_{i}_m_e"] for i in range(5)])

    @property
    def t_3_mo(self):
        return self._t_n_o(*[self.params[f"u_3_{i}_m_o"] for i in range(3)])

    @property
    def t_4_me(self):
        return self._t_n_e(*[self.params[f"u_4_{i}_m_e"] for i in range(5)])

    @property
    def t_4_mo(self):
        return self._t_n_o(*[self.params[f"u_4_{i}_m_o"] for i in range(3)])

    @property
    def t_5_me(self):
        return np.array([
            [self.params["u_5_0_m_e"],  -self.params["u_5_1_m_e"],  0.0],
            [-self.params["u_5_6_m_e"], self.params["u_5_3_m_e"],   0.0],
            [0.0,                       0.0,                        self.params["u_5_5_m_e"]]
        ])

    @property
    def t_5_mo(self):
        return np.diag((self.params["u_5_2_m_o"], self.params["u_5_0_m_o"]))

    @property
    def t_5_xe(self):
        return self._t_5_xr(*[self.params[f"u_5_{i}_x_e"] for i in ("0", "2", "3", "5", "6")])

    @property
    def t_5_xo(self):
        return self._t_5_xr(*[self.params[f"u_5_{i}_x_o"] for i in ("0", "2", "3", "5", "6")])

    @property
    def t_6_me(self):
        return self._t_m_me(*[self.params[f"u_6_{i}_m_e"] for i in range(6)])

    @property
    def t_6_mo(self):
        return self._t_m_mo(*[self.params[f"u_6_{i}_m_o"] for i in range(3)])

    @property
    def t_6_xe(self):
        return self._t_m_xr(*[self.params[f"u_6_{i}_x_e"] for i in range(6)])

    @property
    def t_6_xo(self):
        return self._t_m_xr(*[self.params[f"u_6_{i}_x_o"] for i in range(6)])
