"""Tight-binding models transition metal dichalcogenides"""
from .tmd_abstract_lattice import AbstractLattice, LatticeOrbitals
from .tmd_hopping_matrices import TmdMatrices
from .parameters import fang, liu2, liu6, wu, jorissen, all, cappelluti, dias


class TmdNN2Me(AbstractLattice):
    """3 bands model for even metal orbitals with the second nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2]},
            orbs={"M": ["dz2", "dx2y2", "dxy"]},
            group={"M": [0, 1, 1]}
        )
        super().__init__(orbital=orbital, params=liu2["MoS2"], lattice_name="3 bands 2NN model",
                         n_v=0, n_b=3)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = t_m.e_me
        h_2_m = t_m.t_2_me
        keys = ["h_0_m", "h_2_m", "a", "lamb_m"]
        values = [h_0_m, h_2_m, self.params["a"], self.params["lamb_m"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN12MeXe(AbstractLattice):
    """6 bands model for even metal and chalcogen orbitals with the first- and second-nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2], "X": [1, -1, 0]},
            orbs={"M": ["dz2", "dx2y2", "dxy"], "X": ["pxe", "pye", "pze"]},
            group={"M": [0, 2, 2], "X": [0, 0, 1]}
        )
        super().__init__(orbital=orbital, params=jorissen["MoS2"], lattice_name="6 bands 2NN model",
                         n_v=3, n_b=6)
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
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN12MeoXeo(AbstractLattice):
    """3 bands model for metal and chalocogen orbitals with the first- and second-nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
            orbs={"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"], "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
            group={"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}
        )
        super().__init__(orbital=orbital, params=cappelluti["MoS2"], lattice_name="11 bands 2NN model",
                         n_v=6, n_b=11)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]
    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_0_x = self.block_diag(t_m.e_xe, t_m.e_xo)
        h_1_m = self.block_diag(t_m.t_1_me, t_m.t_1_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_2_x = self.block_diag(t_m.t_2_xe, t_m.t_2_xo)
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN123MeoXeo(AbstractLattice):
    """11 bands model for metal and chalcogen orbitals with the first-, second- and third-nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
            orbs={"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"], "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
            group={"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}
        )
        super().__init__(orbital=orbital, params=fang["MoS2"], lattice_name="11 bands 3NN model",
                         n_v=6, n_b=11)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_0_x = self.block_diag(t_m.e_xe, t_m.e_xo)
        h_1_m = self.block_diag(t_m.t_1_me, t_m.t_1_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_2_x = self.block_diag(t_m.t_2_xe, t_m.t_2_xo)
        h_3_m = self.block_diag(t_m.t_3_me, t_m.t_1_mo * 0)
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "h_3_m",
                "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, h_3_m,
                  self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN125MeoXeo(AbstractLattice):
    """11 bands model for metal and chalcogen orbitals with the first-, second- and fifth- nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
            orbs={"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"], "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
            group={"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}
        )
        super().__init__(orbital=orbital, params=dias["MoS2"], lattice_name="11 bands 5NN model",
                         n_v=6, n_b=11)
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
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "h_5_m", "h_5_c",
                "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, h_5_m, h_5_x,
                  self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN256Me(AbstractLattice):
    """3 bands model for even metal orbitals with the second-, fifth- and sixth-nearest-neighbor hopping"""

    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2]},
            orbs={"M": ["dz2", "dx2y2", "dxy"]},
            group={"M": [0, 1, 1]}
        )
        super().__init__(orbital=orbital, params=liu6["MoS2"], lattice_name="3 bands 6NN model",
                         n_v=0, n_b=3)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]
    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = t_m.e_me
        h_2_m = t_m.t_2_me
        h_5_m = t_m.t_5_me
        h_6_m = t_m.t_6_me
        keys = ["h_0_m", "h_2_m", "h_5_m", "h_6_m", "a", "lamb_m"]
        values = [h_0_m, h_2_m, h_5_m, h_6_m, self.params["a"], self.params["lamb_m"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN256Meo(AbstractLattice):
    """5 bands model for metal orbitals with the second-, fifth- and sixth-nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2, 1, -1]},
            orbs={"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"]},
            group={"M": [0, 1, 1, 2, 2]}
        )
        super().__init__(orbital=orbital, params=wu["MoS2"], lattice_name="5 bands 6NN model",
                         n_v=0, n_b=5)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_5_m = self.block_diag(t_m.t_5_me, t_m.t_5_mo)
        h_6_m = self.block_diag(t_m.t_6_me, t_m.t_6_mo)
        keys = ["h_0_m", "h_2_m", "h_5_m", "h_6_m", "a", "lamb_m"]
        values = [h_0_m, h_2_m, h_5_m, h_6_m, self.params["a"], self.params["lamb_m"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))


class TmdNN123456MeoXeo(AbstractLattice):
    """11 bands model for metal and chalcogen orbitals with all the hoppings up to sixth nearest-neighbor hopping"""
    def __init__(self, **kwargs):
        orbital = LatticeOrbitals(
            l_number={"M": [0, 2, -2, 1, -1], "X": [1, -1, 0, 1, -1, 0]},
            orbs={"M": ["dz2", "dx2y2", "dxy", "dxz", "dyz"], "X": ["pxe", "pye", "pze", "pxo", "pyo", "pzo"]},
            group={"M": [0, 1, 1, 2, 2], "X": [0, 0, 1, 2, 2, 3]}
        )
        super().__init__(orbital=orbital, params=all["MoS2"], lattice_name="11 bands 6NN model",
                         n_v=6, n_b=11)
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    def _generate_matrices(self):
        t_m = TmdMatrices(self.params)
        h_0_m = self.block_diag(t_m.e_me, t_m.e_mo)
        h_0_x = self.block_diag(t_m.e_xe, t_m.e_xo)
        h_1_m = self.block_diag(t_m.t_1_me, t_m.t_1_mo)
        h_2_m = self.block_diag(t_m.t_2_me, t_m.t_2_mo)
        h_2_x = self.block_diag(t_m.t_2_xe, t_m.t_2_xo)
        h_3_m = self.block_diag(t_m.t_3_me, t_m.t_3_mo)
        h_4_m = self.block_diag(t_m.t_4_me, t_m.t_4_mo)
        h_5_m = self.block_diag(t_m.t_5_me, t_m.t_5_mo)
        h_5_x = self.block_diag(t_m.t_5_xe, t_m.t_5_xo)
        h_6_m = self.block_diag(t_m.t_6_me, t_m.t_6_mo)
        h_6_x = self.block_diag(t_m.t_6_xe, t_m.t_6_xo)
        keys = ["h_0_m", "h_0_c", "h_1_m", "h_2_m", "h_2_c", "h_3_m", "h_4_m", "h_5_m", "h_5_c", "h_6_m", "h_6_c",
                "a", "lamb_m", "lamb_c"]
        values = [h_0_m, h_0_x, h_1_m, h_2_m, h_2_x, h_3_m, h_4_m, h_5_m, h_5_x, h_6_m, h_6_x,
                  self.params["a"], self.params["lamb_m"], self.params["lamb_x"]]
        self.lattice_params.set_params(dict([(key, value) for key, value in zip(keys, values)]))
