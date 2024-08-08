"""The matrices for the strained TMD hoppings."""
import numpy as np
from .parameters.strain import StrainParametersList
from .tmd_hopping_matrices import TmdMatrices


class TmdStrainMatrices(TmdMatrices):
    """Construct the strained TMD hopping matrices"""

    def __init__(self, params: StrainParametersList):
        """Initialize the strained TMD hopping matrices

        Parameters:
            params (StrainParametersList): The parameters for the strained TMD lattice."""
        super().__init__(params)
        self.params: StrainParametersList = params

    @staticmethod
    def _e_xr_u(u_0, u_1) -> np.ndarray:
        """The onsite matrix for the chalcogen orbitals under uniaxial strain"""
        return np.array([[u_0, 0, 0], [0, -u_0, u_1], [0, u_1, 0]])

    @staticmethod
    def _e_me_u(u_0, u_1) -> np.ndarray:
        """The onsite matrix for the even metal orbitals under uniaxial strain"""
        return np.array([[0, u_0, 0], [u_0, -u_1, 0], [0, 0, u_1]])

    @staticmethod
    def _e_mo_u(u_0) -> np.ndarray:
        """The onsite matrix for the odd metal orbitals under uniaxial strain"""
        return np.array([[u_0, 0], [0, -u_0]])

    @staticmethod
    def _e_xr_s(u_0, u_1) -> np.ndarray:
        """The onsite matrix for the chalcogen orbitals under shear strain"""
        return np.array([[0, -u_0, -u_1], [-u_0, 0, 0], [-u_1, 0, 0]])

    @staticmethod
    def _e_me_s(u_0, u_1) -> np.ndarray:
        """The onsite matrix for the even metal orbitals under shear strain"""
        return np.array([[0, 0, -u_0], [0, 0, -u_1], [-u_0, -u_1, 0]])

    @staticmethod
    def _e_mo_s(u_0) -> np.ndarray:
        """The onsite matrix for the odd metal orbitals under shear strain"""
        return np.array([[0, -u_0], [-u_0, 0]])

    @staticmethod
    def _t_n_e_s(u_5, u_6, u_7, u_8) -> np.ndarray:
        """The hopping matrices for the first-, third- and fourth-nearest neighbours from the even metal orbitals under shear strain"""
        return np.array([[u_5, u_6, 0], [0, 0, u_7], [0, 0, u_8]])

    @staticmethod
    def _t_n_o_s(u_3, u_4, u_5) -> np.ndarray:
        """The hopping matrices for the first-, third- and fourth-nearest neighbours from the odd metal orbitals under shear strain"""
        return np.array([[0, u_3], [u_4, 0], [u_5, 0]])

    @staticmethod
    def _t_m_xr_s(u_1, u_2, u_4) -> np.ndarray:
        """The hopping matrices for the second- and sixth-nearest neighbours from the chalcogen orbitals under shear strain"""
        return np.array([[0, u_1, u_2], [u_1, 0, u_4], [u_2, -u_4, 0]])

    @staticmethod
    def _t_m_me_s(u_1, u_2, u_4) -> np.ndarray:
        """The hopping matrices for the second-, fifth- and sixth-nearest neighbours from the even metal orbitals under shear strain"""
        return np.array([[0, u_1, u_2], [-u_1, 0, u_4], [u_2, u_4, 0]])

    @staticmethod
    def _t_m_mo_s(u_1) -> np.ndarray:
        """The hopping matrices for the second-, fifth- and sixth-nearest neighbours from the odd metal orbitals under shear strain"""
        return np.array([[0, u_1], [u_1, 0]])

    @staticmethod
    def _t_5_me_s(u_2, u_4, u_7, u_8) -> np.ndarray:
        """The hopping matrix for the fifth-nearest neighbours from the even metal orbitals under shear strain"""
        return np.array([[0, 0, u_2], [0, 0, u_4], [u_7, u_8, 0]])

    @staticmethod
    def _t_5_mo_s(u_1, u_3) -> np.ndarray:
        """The hopping matrix for the fifth-nearest neighbours from the odd metal orbitals under shear strain"""
        return np.array([[0, u_1], [u_3, 0]])

    @staticmethod
    def _t_5_xr_s(u_1, u_4, u_7, u_8) -> np.ndarray:
        """The hopping matrices for the fifth-nearest neighbours from the chalcogen orbitals under shear strain"""
        return np.array([[0, u_1, u_4], [u_7, 0, 0], [u_8, 0, 0]])

    @property
    def e_xe_b(self) -> list[np.ndarray]:
        """The hopping matrix for the even chalcogen orbitals under biaxial strain"""
        return [
            np.diag((
                self.params[f"eps_0_x_e_b_{order_idx}"],
                self.params[f"eps_0_x_e_b_{order_idx}"],
                self.params[f"eps_1_x_e_b_{order_idx}"]
            ))
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_xo_b(self) -> list[np.ndarray]:
        """The hopping matrix for the odd chalcogen orbitals under biaxial strain"""
        return [
            np.diag((
                self.params[f"eps_0_x_o_b_{order_idx}"],
                self.params[f"eps_0_x_o_b_{order_idx}"],
                self.params[f"eps_1_x_o_b_{order_idx}"]
            ))
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_me_b(self) -> list[np.ndarray]:
        """The hopping matrix for the even metal orbitals under biaxial strain"""
        return [
            np.diag((
                self.params[f"eps_0_m_e_b_{order_idx}"],
                self.params[f"eps_1_m_e_b_{order_idx}"],
                self.params[f"eps_1_m_e_b_{order_idx}"]
            ))
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_mo_b(self) -> list[np.ndarray]:
        """The hopping matrix for the odd metal orbitals under biaxial strain"""
        return [
            np.diag((
                self.params[f"eps_0_m_e_b_{order_idx}"],
                self.params[f"eps_0_m_e_b_{order_idx}"]
            ))
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_xe_u(self) -> list[np.ndarray]:
        """The hopping matrix for the even chalcogen orbitals under uniaxial strain"""
        return [self._e_xr_u(
            self.params[f"eps_0_x_e_u_{order_idx}"],
            self.params[f"eps_1_x_e_u_{order_idx}"]
        ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_xo_u(self) -> list[np.ndarray]:
        """The hopping matrix for the odd chalcogen orbitals under uniaxial strain"""
        return [
            self._e_xr_u(
                self.params[f"eps_0_x_o_u_{order_idx}"],
                self.params[f"eps_1_x_o_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_me_u(self) -> list[np.ndarray]:
        """The hopping matrix for the even metal orbitals under uniaxial strain"""
        return [
            self._e_me_u(
                self.params[f"eps_0_m_e_u_{order_idx}"],
                self.params[f"eps_1_m_e_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_mo_u(self) -> list[np.ndarray]:
        """The hopping matrix for the odd metal orbitals under uniaxial strain"""
        return [
            self._e_mo_u(
                self.params[f"eps_0_m_o_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_xe_s(self) -> list[np.ndarray]:
        """The hopping matrix for the even chalcogen orbitals under shear strain"""
        return [
            self._e_xr_s(
                self.params[f"eps_0_x_e_u_{order_idx}"],
                self.params[f"eps_1_x_e_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_xo_s(self) -> list[np.ndarray]:
        """The hopping matrix for the odd chalcogen orbitals under shear strain"""
        return [
            self._e_xr_s(
                self.params[f"eps_0_x_o_u_{order_idx}"],
                self.params[f"eps_1_x_o_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_me_s(self) -> list[np.ndarray]:
        """The hopping matrix for the even metal orbitals under shear strain"""
        return [
            self._e_me_s(
                self.params[f"eps_0_m_e_u_{order_idx}"],
                self.params[f"eps_1_m_e_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def e_mo_s(self) -> list[np.ndarray]:
        """The hopping matrix for the odd metal orbitals under shear strain"""
        return [
            self._e_mo_s(
                self.params[f"eps_0_m_o_u_{order_idx}"]
            )
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_me_b(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_1_{i}_m_e_b_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_mo_b(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_1_{i}_m_o_b_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_me_b(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_m_me(*[self.params[f"u_2_{i}_m_e_b_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_mo_b(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_m_mo(*[self.params[f"u_2_{i}_m_o_b_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xe_b(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even chalcogen orbitals under biaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_2_{i}_x_e_b_{order_idx}"] for i in range(6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xo_b(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd chalcogen orbitals under biaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_2_{i}_x_o_b_{order_idx}"] for i in range(6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_me_b(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_3_{i}_m_e_b_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_mo_b(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_3_{i}_m_o_b_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_me_b(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_4_{i}_m_e_b_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_mo_b(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_4_{i}_m_o_b_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_me_b(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_5_me(*[self.params[f"u_5_{i}_m_e_b_{order_idx}"] for i in ("0", "1", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_mo_b(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_5_mo(*[self.params[f"u_5_{i}_m_o_b_{order_idx}"] for i in ("0", "2")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xe_b(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even chalcogen orbitals under biaxial strain"""
        return [
            self._t_5_xr(*[self.params[f"u_5_{i}_x_e_b_{order_idx}"] for i in ("0", "2", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xo_b(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd chalcogen orbitals under biaxial strain"""
        return [
            self._t_5_xr(*[self.params[f"u_5_{i}_x_o_b_{order_idx}"] for i in ("0", "2", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_me_b(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even metal orbitals under biaxial strain"""
        return [
            self._t_m_me(*[self.params[f"u_6_{i}_m_e_b_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_mo_b(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd metal orbitals under biaxial strain"""
        return [
            self._t_m_mo(*[self.params[f"u_6_{i}_m_o_b_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xe_b(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even chalcogen orbitals under biaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_6_{i}_x_e_b_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xo_b(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd chalcogen orbitals under biaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_6_{i}_x_o_b_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_me_u(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_1_{i}_m_e_u_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_mo_u(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_1_{i}_m_o_u_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_me_u(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_m_me(*[self.params[f"u_2_{i}_m_e_u_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_mo_u(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_m_mo(*[self.params[f"u_2_{i}_m_o_u_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xe_u(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even chalcogen orbitals under uniaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_2_{i}_x_e_u_{order_idx}"] for i in range(6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xo_u(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd chalcogen orbitals under uniaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_2_{i}_x_o_u_{order_idx}"] for i in range(6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_me_u(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_3_{i}_m_e_u_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_mo_u(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_3_{i}_m_o_u_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_me_u(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_n_e(*[self.params[f"u_4_{i}_m_e_u_{order_idx}"] for i in range(5)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_mo_u(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_n_o(*[self.params[f"u_4_{i}_m_o_u_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_me_u(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_5_me(*[self.params[f"u_5_{i}_m_e_u_{order_idx}"] for i in ("0", "1", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_mo_u(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_5_mo(*[self.params[f"u_5_{i}_m_o_u_{order_idx}"] for i in ("0", "2")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xe_u(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even chalcogen orbitals under uniaxial strain"""
        return [
            self._t_5_xr(*[self.params[f"u_5_{i}_x_e_u_{order_idx}"] for i in ("0", "2", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xo_u(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd chalcogen orbitals under uniaxial strain"""
        return [
            self._t_5_xr(*[self.params[f"u_5_{i}_x_o_u_{order_idx}"] for i in ("0", "2", "3", "5", "6")])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_me_u(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even metal orbitals under uniaxial strain"""
        return [
            self._t_m_me(*[self.params[f"u_6_{i}_m_e_u_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_mo_u(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd metal orbitals under uniaxial strain"""
        return [
            self._t_m_mo(*[self.params[f"u_6_{i}_m_o_u_{order_idx}"] for i in range(3)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xe_u(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even chalcogen orbitals under uniaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_6_{i}_x_e_u_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xo_u(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd chalcogen orbitals under uniaxial strain"""
        return [
            self._t_m_xr(*[self.params[f"u_6_{i}_x_o_u_{order_idx}"] for i in range(6)])
            for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_me_s(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_n_e_s(*[self.params[f"u_1_{i}_m_e_s_{order_idx}"] for i in range(5, 9)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_1_mo_s(self) -> list[np.ndarray]:
        """The first-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_n_o_s(*[self.params[f"u_1_{i}_m_o_s_{order_idx}"] for i in range(3, 6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_me_s(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_m_me_s(*[self.params[f"u_2_{i}_m_e_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_mo_s(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_m_mo_s(*[self.params[f"u_2_{i}_m_o_s_{order_idx}"] for i in [1]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xe_s(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the even chalcogen orbitals under shear strain"""
        return [
            self._t_m_xr_s(*[self.params[f"u_2_{i}_x_e_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_2_xo_s(self) -> list[np.ndarray]:
        """The second-nearest neighbour hopping matrix from the odd chalcogen orbitals under shear strain"""
        return [
            self._t_m_xr_s(*[self.params[f"u_2_{i}_x_o_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_me_s(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_n_e_s(*[self.params[f"u_3_{i}_m_e_s_{order_idx}"] for i in range(5, 9)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_3_mo_s(self) -> list[np.ndarray]:
        """The third-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_n_o_s(*[self.params[f"u_3_{i}_m_o_s_{order_idx}"] for i in range(3, 6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_me_s(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_n_e_s(*[self.params[f"u_4_{i}_m_e_s_{order_idx}"] for i in range(5, 9)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_4_mo_s(self) -> list[np.ndarray]:
        """The fourth-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_n_o_s(*[self.params[f"u_4_{i}_m_o_s_{order_idx}"] for i in range(3, 6)]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_me_s(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_5_me_s(*[self.params[f"u_5_{i}_m_e_s_{order_idx}"] for i in [2, 4, 7, 8]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_mo_s(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_5_mo_s(*[self.params[f"u_5_{i}_m_o_s_{order_idx}"] for i in [1, 3]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xe_s(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the even chalcogen orbitals under shear strain"""
        return [
            self._t_5_xr_s(*[self.params[f"u_5_{i}_x_e_s_{order_idx}"] for i in [1, 4, 7, 8]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_5_xo_s(self) -> list[np.ndarray]:
        """The fifth-nearest neighbour hopping matrix from the odd chalcogen orbitals under shear strain"""
        return [
            self._t_5_xr_s(*[self.params[f"u_5_{i}_x_o_s_{order_idx}"] for i in [1, 4, 7, 8]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_me_s(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even metal orbitals under shear strain"""
        return [
            self._t_m_me_s(*[self.params[f"u_6_{i}_m_e_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_mo_s(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd metal orbitals under shear strain"""
        return [
            self._t_m_mo_s(*[self.params[f"u_6_{i}_m_o_s_{order_idx}"] for i in [1]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xe_s(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the even chalcogen orbitals under shear strain"""
        return [
            self._t_m_xr_s(*[self.params[f"u_6_{i}_x_e_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]

    @property
    def t_6_xo_s(self) -> list[np.ndarray]:
        """The sixth-nearest neighbour hopping matrix from the odd chalcogen orbitals under shear strain"""
        return [
            self._t_m_xr_s(*[self.params[f"u_6_{i}_x_o_s_{order_idx}"] for i in [1, 2, 4]]
            ) for order_idx in range(1, self.params.max_order+1)]
