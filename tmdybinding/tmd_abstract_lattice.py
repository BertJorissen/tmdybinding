"""Tight-binding models for group 1 transition metal dichalcogenides (tmd), 6 band."""
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import re
from typing import List, Optional, Tuple
from .parameters import ParametersList


class LatticeParams:
    """
    Class for saving the matrices in the abstract-tmd-lattice-contructor-class. There are also checks for a consistent
    model.
    Usage: construct an instance, and call this instance or the set_params function, with the arguments and values
    the matrices and values of these respectively. If the set_params function is called again, the params are set to
    None if not gooed defined. The argument "keep=True" can be passed to keep the old arguments.
    """
    def __init__(self, **kwargs):
        self.__h_0_m = None
        self.__h_0_c = None
        self.__h_1_m = None
        self.__h_2_m = None
        self.__h_2_c = None
        self.__h_3_m = None
        self.__h_4_m = None
        self.__h_5_m = None
        self.__h_5_c = None
        self.__h_6_m = None
        self.__h_6_c = None
        self.__a = None
        self.__lamb_m = None
        self.__lamb_c = None
        self.__keys = [
            "h_0_m", "h_0_c",
            "h_1_m",
            "h_2_m", "h_2_c",
            "h_3_m",
            "h_4_m",
            "h_5_m", "h_5_c",
            "h_6_m", "h_6_c",
            "a", "lamb_m", "lamb_c"
        ]
        if not kwargs == {}:
            self.set_params(**kwargs)

    def __call__(self, **kwargs):
        self.set_params(**kwargs)

    @staticmethod
    def _params_check(params, h_x):
        return (np.array(params[h_x]) if params[h_x] is not None else None) if h_x in [*params] else None

    def set_params(self, **kwargs):
        params = self._keys_dict if "keep" in [*kwargs] and kwargs.pop("keep") else {}
        for name in [*kwargs]:
            assert name in self.__keys, "key not in expected hoppings, possible wrong name"
            params[name] = kwargs[name]
        self._component_check(params)
        self._shape_check(params)
        self.__h_0_m = self._params_check(params, "h_0_m")
        self.__h_0_c = self._params_check(params, "h_0_c")
        self.__h_1_m = self._params_check(params, "h_1_m")
        self.__h_2_m = self._params_check(params, "h_2_m")
        self.__h_2_c = self._params_check(params, "h_2_c")
        self.__h_3_m = self._params_check(params, "h_3_m")
        self.__h_4_m = self._params_check(params, "h_4_m")
        self.__h_5_m = self._params_check(params, "h_5_m")
        self.__h_5_c = self._params_check(params, "h_5_c")
        self.__h_6_m = self._params_check(params, "h_6_m")
        self.__h_6_c = self._params_check(params, "h_6_c")
        self.__a = self._params_check(params, "a")
        self.__lamb_m = self._params_check(params, "lamb_m")
        self.__lamb_c = self._params_check(params, "lamb_c")

    def _make_bools(self, kwargs):
        m_bool = self._attr_check("h_0_m", kwargs)
        c_bool = self._attr_check("h_0_c", kwargs)

        h_1_m_bool = self._attr_check("h_1_m", kwargs)  # M-X
        h_2_m_bool = self._attr_check("h_2_m", kwargs)  # M-M
        h_2_c_bool = self._attr_check("h_2_c", kwargs)  # X-X
        h_3_m_bool = self._attr_check("h_3_m", kwargs)  # M-X
        h_4_m_bool = self._attr_check("h_4_m", kwargs)  # M-X
        h_5_m_bool = self._attr_check("h_5_m", kwargs)  # M-M
        h_5_c_bool = self._attr_check("h_5_c", kwargs)  # X-X
        h_6_m_bool = self._attr_check("h_6_m", kwargs)  # M-M
        h_6_c_bool = self._attr_check("h_6_c", kwargs)  # M-M
        return (m_bool, c_bool,
                h_1_m_bool,
                h_2_m_bool, h_2_c_bool,
                h_3_m_bool,
                h_4_m_bool,
                h_5_m_bool, h_5_c_bool,
                h_6_m_bool, h_6_c_bool)

    def _component_check(self, kwargs):
        (m_bool, c_bool,
         h_1_m_bool,
         h_2_m_bool, h_2_c_bool,
         h_3_m_bool,
         h_4_m_bool,
         h_5_m_bool, h_5_c_bool,
         h_6_m_bool, h_6_c_bool) = self._make_bools(kwargs)
        assert m_bool or c_bool, "not X nor M in model"
        assert not (not m_bool and h_1_m_bool), "not M in model and 1st hopping from M"
        assert not (not c_bool and h_1_m_bool), "not M in model and 1st hopping to X"
        assert not (not m_bool and h_2_m_bool), "not M in model and 2nd hopping from/to M"
        assert not (not c_bool and h_2_c_bool), "not M in model and 2nd hopping from/to X"
        assert not (not m_bool and h_3_m_bool), "not M in model and 3rd hopping from M"
        assert not (not c_bool and h_3_m_bool), "not M in model and 3rd hopping to X"
        assert not (not m_bool and h_4_m_bool), "not M in model and 4rd hopping from M"
        assert not (not c_bool and h_4_m_bool), "not M in model and 4rd hopping to X"
        assert not (not m_bool and h_5_m_bool), "not M in model and 5th hopping from/to M"
        assert not (not c_bool and h_5_c_bool), "not M in model and 5th hopping from/to X"
        assert not (not m_bool and h_6_m_bool), "not M in model and 6th hopping from/to M"
        assert not (not c_bool and h_6_c_bool), "not M in model and 6th hopping from/to X"
        assert h_1_m_bool or h_2_c_bool or h_2_m_bool or h_3_m_bool or h_4_m_bool or h_5_c_bool or h_5_m_bool or\
               h_6_c_bool or h_6_m_bool, "no hoppings specified in the model"

    def _shape_check(self, kwargs):
        (m_bool, c_bool,
         h_1_m_bool,
         h_2_m_bool, h_2_c_bool,
         h_3_m_bool,
         h_4_m_bool,
         h_5_m_bool, h_5_c_bool,
         h_6_m_bool, h_6_c_bool) = self._make_bools(kwargs)
        if m_bool:
            m_shape = self._check_shape(kwargs["h_0_m"])
            if h_1_m_bool:
                assert np.shape(kwargs["h_1_m"])[1] == m_shape, "shape 1st hopping from M not correct"
            if h_2_m_bool:
                assert np.shape(kwargs["h_2_m"])[0] == m_shape, "shape 2nd hopping to M not correct"
                assert np.shape(kwargs["h_2_m"])[1] == m_shape, "shape 2nd hopping from M not correct"
            if h_3_m_bool:
                assert np.shape(kwargs["h_3_m"])[1] == m_shape, "shape 3rd hopping from M not correct"
            if h_4_m_bool:
                assert np.shape(kwargs["h_4_m"])[1] == m_shape, "shape 4rd hopping from M not correct"
            if h_5_m_bool:
                assert np.shape(kwargs["h_5_m"])[0] == m_shape, "shape 5th hopping to M not correct"
                assert np.shape(kwargs["h_5_m"])[1] == m_shape, "shape 5th hopping from M not correct"
            if h_6_m_bool:
                assert np.shape(kwargs["h_6_m"])[0] == m_shape, "shape 6th hopping to M not correct"
                assert np.shape(kwargs["h_6_m"])[1] == m_shape, "shape 6th hopping from M not correct"
        if c_bool:
            c_shape = self._check_shape(kwargs["h_0_c"])
            if h_1_m_bool:
                assert np.shape(kwargs["h_1_m"])[0] == c_shape, "shape 1st hopping to X not correct"
            if h_2_c_bool:
                assert np.shape(kwargs["h_2_c"])[0] == c_shape, "shape 2nd hopping to X not correct"
                assert np.shape(kwargs["h_2_c"])[1] == c_shape
            if h_3_m_bool:
                assert np.shape(kwargs["h_3_m"])[0] == c_shape, "shape 3rd hopping to X not correct"
            if h_4_m_bool:
                assert np.shape(kwargs["h_4_m"])[0] == c_shape, "shape 4rd hopping to X not correct"
            if h_5_c_bool:
                assert np.shape(kwargs["h_5_c"])[0] == c_shape, "shape 5th hopping to X not correct"
                assert np.shape(kwargs["h_5_c"])[1] == c_shape, "shape 5th hopping from X not correct"
            if h_6_c_bool:
                assert np.shape(kwargs["h_6_c"])[0] == c_shape, "shape 6th hopping to X not correct"
                assert np.shape(kwargs["h_6_c"])[1] == c_shape, "shape 6th hopping from X not correct"

    @staticmethod
    def _attr_check(name, params):
        if name in params:
            return params[name] is not None
        else:
            return False

    @staticmethod
    def _check_shape(h_0):
        h_0_shape = np.shape(h_0)
        assert len(h_0_shape) < 3, "onsite energy too many dims (max 2 as matrix)"
        if len(h_0_shape) == 2:
            assert h_0_shape[0] == h_0_shape[1] or h_0_shape[0] == 1, "shape of x and y not the same for onsite energy"
            return h_0_shape[1]
        else:
            return h_0_shape[0]

    @property
    def _keys_dict(self):
        dict_list = []
        for name in self.__keys:
            value = getattr(self, name)
            if value is not None:
                dict_list.append((name, value))
        return dict(dict_list)

    @property
    def h_0_m(self):
        return self.__h_0_m

    @property
    def h_0_c(self):
        return self.__h_0_c

    @property
    def h_1_m(self):
        return self.__h_1_m

    @property
    def h_2_m(self):
        return self.__h_2_m

    @property
    def h_2_c(self):
        return self.__h_2_c

    @property
    def h_3_m(self):
        return self.__h_3_m

    @property
    def h_4_m(self):
        return self.__h_4_m

    @property
    def h_5_m(self):
        return self.__h_5_m

    @property
    def h_5_c(self):
        return self.__h_5_c

    @property
    def h_6_m(self):
        return self.__h_6_m

    @property
    def h_6_c(self):
        return self.__h_6_c

    @property
    def a(self):
        return self.__a

    @property
    def lamb_m(self):
        return self.__lamb_m

    @property
    def lamb_c(self):
        return self.__lamb_c


class LatticeOrbitals:
    """
    object for saving the l and s values of a lattice model. The input are the relevent parameters given in a dict
    with the right name of the orbital.
    """
    def __init__(self, **kwargs):
        self.__l = None
        self.__orbs = None
        self.__names = None
        self.__group = None
        self.clockwise = kwargs.pop("clockwise") if "clockwise" in [*kwargs] else False
        if not kwargs == {}:
            self.set_params(**kwargs)

    def set_params(self, **kwargs):
        self._check_names(kwargs)
        self._check_shape(kwargs)

    @staticmethod
    def _attr_check(name, params):
        if name in params:
            return params[name] is not None
        else:
            return False

    def _make_bools(self, kwargs):
        l_bool = self._attr_check("l", kwargs)
        orbs_bool = self._attr_check("orbs", kwargs)
        group_bool = self._attr_check("group", kwargs)
        assert l_bool or orbs_bool, "not enough input given, at least l or orbs needed"
        return l_bool, orbs_bool, group_bool

    def _check_names(self, kwargs):
        (l_bool, orbs_bool, group_bool) = self._make_bools(kwargs)
        names = None
        if l_bool:
            if names is None:
                names = [*kwargs['l']]
        if orbs_bool:
            if names is None:
                names = [*kwargs['orbs']]
            else:
                assert names == [*kwargs['orbs']], "the names of orbs and l are not the same"
        if group_bool:
            if names is None:
                names = [*kwargs['group']]
            else:
                assert names == [*kwargs['group']], "the names of group and l and/or orbs are not the same"
        self.__names = names

    def _check_shape(self, kwargs):
        # see if shape is similar
        (l_bool, orbs_bool, group_bool) = self._make_bools(kwargs)
        shape = None
        if l_bool:
            if shape is None:
                shape = [np.shape(kwargs['l'][name]) for name in self.names]
        if orbs_bool:
            if shape is None:
                shape = [np.shape(kwargs['orbs'][name]) for name in self.names]
            else:
                assert shape == [np.shape(kwargs['orbs'][name]) for name in self.names],\
                    "the shape of orbs and l are not the same"
        if group_bool:
            if shape is None:
                shape = [np.shape(kwargs['group'][name]) for name in self.names]
            else:
                assert shape == [np.shape(kwargs['group'][name]) for name in self.names],\
                    "the shape of group and l and/or orbs are not the same"
        # see if definition of group and l are correct. If no group is specified, but there is an l, the l is checked
        if not group_bool and l_bool:
            for name in self.names:
                l = np.array(kwargs["l"][name])
                for lm in set(np.abs(l)):
                    l_b = np.abs(l) == lm
                    assert np.sum(l_b) < 3, "the representation can only be given for max two parts"
                    assert (np.sum(l_b) == 2 and l[l_b][0] == -l[l_b][1]) or np.sum(l_b) == 1,\
                        "can't have same l-number if group is not defined"
        # create l or group is one of them is missing
        l_out = kwargs['l'] if l_bool else dict([(name, np.zeros(shape[j])) for j, name in enumerate(self.names)])
        group_out = kwargs['group'] if group_bool else (
             dict([(name, np.abs(l_out[name])) for name in self.names]) if l_bool else
             dict([(name, np.arange(shape[j])) for j, name in enumerate(self.names)]))
        # check l and group
        for name in self.names:
            group_l = np.array(group_out[name])
            l = np.array(l_out[name])
            for group_i in set(group_l):
                group_b = group_i == group_l
                assert np.sum(group_b) < 3, "the group for representation can only be given for max two parts"
                assert (np.sum(group_b) == 2 and l[group_b][0] == -l[group_b][1]) or np.sum(group_b) == 1, \
                    "can't have same l-number in same group"
        self.__l = l_out
        self.__orbs = kwargs['orbs'] if orbs_bool else dict([(name, [str(k) for k in np.arange(shape[j])])
                                                             for j, name in enumerate(self.names)])
        self.__group = group_out

    @staticmethod
    def rot_mat(phi=0):
        return np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

    @property
    def l(self):
        return self.__l

    @property
    def names(self):
        return self.__names

    @property
    def orbs(self):
        return self.__orbs

    @property
    def group(self):
        return self.__group

    def _make_matrix(self, matrix_func, single):
        dict_list = []
        for name in self.names:
            group_l = np.array(self.group[name])
            l = np.array(self.l[name])
            value = np.zeros((len(l), len(l)))
            for group_i in set(group_l):
                group_b = group_i == group_l
                lm = np.abs(l[group_b][0])
                sign = -1 if l[group_b][0] < 0 else 1
                if np.sum(group_b) == 2:
                    matrix = matrix_func(lm, sign)
                    for ik, k in enumerate(np.arange(len(l))[group_b]):
                        for ikk, kk in enumerate(np.arange(len(l))[group_b]):
                            value[k, kk] = matrix[ik, ikk]
                else:
                    value[np.arange(len(l))[group_b], np.arange(len(l))[group_b]] = single(lm, sign)
            dict_list.append((name, value))
        return dict(dict_list)

    @property
    def ur(self):
        def matrix_func(lm, sign):
            ur_t = self.rot_mat(np.pi * 2 / 3 * lm)
            if sign == -1:
                ur_t = np.transpose(ur_t)
            if self.clockwise:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 1

        return self._make_matrix(matrix_func, single)

    @property
    def sr(self):
        """Mirror on yz-plane"""
        def matrix_func(lm, sign):
            return sign * np.diag([-1, 1]) * (1 if np.abs(lm) == 1 else -1)

        def single(lm, sign):
            return 1

        return self._make_matrix(matrix_func, single)

    def ur_angle(self, angle):
        def matrix_func(lm, sign):
            ur_t = self.rot_mat(angle * lm)
            if sign == -1:
                ur_t = np.transpose(ur_t)
            if self.clockwise:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 1

        return self._make_matrix(matrix_func, single)

    @property
    def s_h(self):
        def matrix_func(lm, sign):
            ur_t = lm / 2 * np.array([[0, -1], [1, 0]])
            if sign == -1:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 0

        return self._make_matrix(matrix_func, single)


class AbstractLattice:
    r"""Abstract class for a tmd-lattice

    Parameters
    ----------
    name : str
        Name of the tmd to model.
    params : Optional[dict]
        Replace or add new material parameters. The dictionary entries must be
        in the format `"name": [ a,   ...]
    sz : float
        Z-component of the spin degree of freedom
    """
    def __init__(self, **kwargs):
        self.lattice_params = LatticeParams()
        self.orbital: LatticeOrbitals = LatticeOrbitals(**kwargs.pop("orbital") if "orbital" in [*kwargs] else None)
        self.__n_valence_band: int = kwargs.pop("n_v") if "n_v" in [*kwargs] else 0
        self.__n_bands: int = kwargs.pop("n_b") if "n_b" in [*kwargs] else 0
        self.__name: str = "MoS2"
        self.__x_name: str = "S"
        self.__x_type: str = "X"
        self.__m_name: str = "Mo"
        self.__m_type: str = "M"
        self.__params: ParametersList = ParametersList()
        self.__soc_eo_flip: bool = False
        self.__lattice_name: str = "Abstract Lattice"
        self._lat4: bool = False
        self._even_odd: bool = False
        self.single_orbital: bool = False
        self.soc: bool = False
        self.soc_polarized: bool = False
        self.soc_sz_part: bool = True
        self.sz: float = 1.
        [setattr(self, var, kwargs[var]) for var in [*kwargs]]

    @property
    def soc_eo_flip(self) -> bool:
        return self.__soc_eo_flip

    @soc_eo_flip.setter
    def soc_eo_flip(self, pol_bool: bool):
        if pol_bool:
            if self.lattice_params.h_0_m is not None:
                assert len(self.orbital.l[self.orb_type(self.m_name)]) == 5,\
                    "the metal doesn't have the right shape for spin-flip term"
                assert sorted(self.orbital.l[self.orb_type(self.m_name)]) == sorted([0, 2, -2, 1, -1]),\
                    "the metal l is wrong for spin-flip"
            if self.lattice_params.h_0_c is not None:
                assert len(self.orbital.l[self.orb_type(self.x_name)]) == 6,\
                    "the chal. doesn't have the right shape for spin-flip term"
                assert sorted(self.orbital.l[self.orb_type(self.x_name)]) == sorted([1, -1, 0, 1, -1, 0]),\
                    "the chal. l is wrong for spin-flip"
        self.__soc_eo_flip = pol_bool

    @property
    def params(self) -> ParametersList:
        return self.__params

    @params.setter
    def params(self, params: ParametersList):
        self.__params = params
        self.name = self.params["material"]
        self._generate_matrices()

    @property
    def lattice_name(self) -> str:
        return self.__lattice_name

    @lattice_name.setter
    def lattice_name(self, lattice_name: str):
        self.__lattice_name = lattice_name

    @property
    def name(self) -> str:
        return self.__name

    @property
    def m_name(self) -> str:
        return self.__m_name

    def orb_type(self, z_name: str) -> Optional[str]:
        if z_name == self.m_name or (self.lat4 and z_name == self.m_name + "2"):
            return self.__m_type
        elif z_name == self.x_name or (self.lat4 and z_name == self.x_name + "2"):
            return self.__x_type
        else:
            return None

    def z_name(self, orb_type: str) -> Optional[str]:
        if orb_type == self.__m_type:
            return self.m_name
        elif orb_type == self.__x_type:
            return self.x_name
        else:
            return None

    @property
    def x_name(self) -> str:
        return self.__x_name

    @name.setter
    def name(self, name: str):
        self.__name = name
        self.__m_name, self.__x_name = re.findall("[A-Z][a-z]*", self.name)
        self._generate_matrices()

    @property
    def lat4(self) -> bool:
        return self._lat4

    @lat4.setter
    def lat4(self, lat4: bool):
        self._lat4 = lat4
        self._generate_matrices()

    @property
    def a1(self) -> np.ndarray:
        return np.array([1, 0]) * self.lattice_params.a

    @property
    def a2(self) -> np.ndarray:
        return (np.array([0, np.sqrt(3)]) if self.lat4 else np.array([-1 / 2, np.sqrt(3) / 2])) * self.lattice_params.a

    @property
    def soc_doubled_ham(self) -> bool:
        return not self.soc_polarized if self.soc else False

    @property
    def soc_eo_flip_used(self) -> bool:
        return self.soc and self.soc_eo_flip and not self.soc_polarized

    @property
    def n_valence_band(self) -> int:
        return (self.__n_valence_band + 1) * 2 - 1 if self.soc_doubled_ham else self.__n_valence_band

    @property
    def n_bands(self) -> int:
        return self.__n_bands * (2 if self.soc_doubled_ham else 1)

    @property
    def _m_orbs(self) -> List[str]:
        m_orbs = [self.m_name + mi for mi in self.orbital.orbs[self.orb_type(self.m_name)]]
        if self.soc_doubled_ham:
            m_orbs = [mi + "u" for mi in m_orbs] + [mi + "d" for mi in m_orbs]
        return m_orbs

    @property
    def _m2_orbs(self) -> List[str]:
        m2_orbs = [self.m_name + "2" + mi for mi in self.orbital.orbs[self.orb_type(self.m_name)]]
        if self.soc_doubled_ham:
            m2_orbs = [mi + "u" for mi in m2_orbs] + [mi + "d" for mi in m2_orbs]
        return m2_orbs

    @property
    def _c_orbs(self) -> List[str]:
        c_orbs = [self.x_name + ci for ci in self.orbital.orbs[self.orb_type(self.x_name)]]
        if self.soc_doubled_ham:
            c_orbs = [ci + "u" for ci in c_orbs] + [ci + "d" for ci in c_orbs]
        return c_orbs

    @property
    def _c2_orbs(self) -> List[str]:
        c2_orbs = [self.x_name + "2" + ci for ci in self.orbital.orbs[self.orb_type(self.x_name)]]
        if self.soc_doubled_ham:
            c2_orbs = [ci + "u" for ci in c2_orbs] + [ci + "d" for ci in c2_orbs]
        return c2_orbs

    @staticmethod
    def _make_name(h_name: str, n_i: int, nfi: str, ntj: str) -> str:
        return h_name + "-" + str(n_i) + "-" + nfi + "-" + ntj

    @staticmethod
    def _separate_name(h_name_i) -> Tuple[str, int, str, str]:
        out = re.findall("[^-]+(?:[^-]*)*", h_name_i)
        assert len(out) == 4, "The given string isn't generated by the right function, the length isn't 4"
        out[1] = int(out[1])
        return tuple(out)

    def _generate_matrices(self):
        a = self.params[self.name]
        keys = ["a"]
        values = [a]
        self.lattice_params(**dict([(key, value) for key, value in zip(keys, values)]))

    @staticmethod
    def _ham(h: np.ndarray, ur_l: np.ndarray, ur_r: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return h, ur_l.T.dot(h.dot(ur_r)), ur_l.dot(h.dot(ur_r.T))

    @staticmethod
    def _reorder(matrix: np.ndarray, keys: Tuple[List[int], List[int]]) -> np.ndarray:
        return np.array([[matrix[xi, yi] for yi in keys[1]] for xi in keys[0]])

    def block_diag(self, *matrices) -> np.ndarray:
        if len(matrices) == 1:
            return matrices[0]
        else:
            matrix0 = matrices[0]
            matrix1 = self.block_diag(*matrices[1:])
            shape_0 = np.shape(matrix0)
            shape_1 = np.shape(matrix1)
            z_1 = [shape_0[-2], shape_1[-1]]
            z_2 = [shape_1[-2], shape_0[-1]]
            nd = len(shape_0)
            assert nd == len(shape_1), "The matrices don't have the right sizes in the additional dimensions"
            if nd > 2:
                for sh in shape_1[:-2]:
                    z_1.insert(0, sh)
                    z_2.insert(0, sh)
            return np.concatenate((np.concatenate((matrix0, np.zeros(z_1)), axis=nd-1),
                                   np.concatenate((np.zeros(z_2), matrix1), axis=nd-1)), axis=nd-2)

    def _make_onsite(self, matrix, name, lamb):
        def ham_sz(sz):
            if self.soc_sz_part:
                s_part = sz * lamb * 1j * self.orbital.s_h[self.orb_type(name)]
                return matrix + s_part
            else:
                return np.array(matrix, dtype=complex)

        if self.soc:
            if self.soc_polarized:
                return ham_sz(self.sz)
            else:
                return self.block_diag(ham_sz(self.sz), ham_sz(-self.sz))
        else:
            return ham_sz(0.)

    def _make_h(self, matrix, from_name, to_name):
        (h_1, h_2, h_3) = self._ham(matrix,
                                    self.orbital.ur[self.orb_type(to_name)],
                                    self.orbital.ur[self.orb_type(from_name)])
        if self.soc_doubled_ham:
            h_1 = np.kron(np.eye(2), h_1)
            h_2 = np.kron(np.eye(2), h_2)
            h_3 = np.kron(np.eye(2), h_3)
        return h_1, h_2, h_3

    def _make_h_angle(self, matrix, from_name, to_name, angle):
        (h_1, h_2, h_3) = self._ham(matrix,
                                    self.orbital.ur_angle(angle)[self.orb_type(to_name)],
                                    self.orbital.ur_angle(angle)[self.orb_type(from_name)])
        return h_1, h_2, h_3

    def _add_hopping(self, lat, mat, f_n, t_n, cos, fnl, tnl, h_name):
        hn = self._make_h(mat, f_n, t_n)
        n_n_n = [0, 1, 2]
        if self.lat4:
            hn = [hn[int(ih/2)] for ih in range(6)]
            n_n_n = [0, 0, 1, 1, 2, 2]
        n_n = 6 if self.lat4 else 3
        if self.single_orbital:
            for f_i, nf_i in enumerate(fnl):
                for t_j, nt_j in enumerate(tnl):
                    h_names = [self._make_name(h_name, n_i, nfi, ntj) for n_i, nfi, ntj in zip(n_n_n, nf_i, nt_j)]
                    lat.register_hopping_energies(dict([(h_n_i, h.T[f_i, t_j]) for h_n_i, h in zip(h_names, hn)]))
                    lat.add_hoppings(*[(co, nfi, ntj, h_n_i) for co, h_n_i, nfi, ntj in zip(cos, h_names, nf_i, nt_j)])
        else:
            h_names = [self._make_name(h_name, n_i, nfi, ntj) for n_i, nfi, ntj in zip(n_n_n, fnl, tnl)]
            lat.register_hopping_energies(dict([(h_n_i, h.T) for h_n_i, h in zip(h_names, hn)]))
            lat.add_hoppings(*[(co, nfi, ntj, h_n_i) for co, h_n_i, nfi, ntj in zip(cos, h_names, fnl, tnl)])
        return lat

    def lattice(self):
        lat = pb.Lattice(a1=self.a1, a2=self.a2)

        if self.lattice_params.h_0_m is not None:
            m_orbs = self._m_orbs
            m2_orbs = self._m2_orbs
            h_0_m = self._make_onsite(self.lattice_params.h_0_m, self.m_name, self.lattice_params.lamb_m)
            if self.soc_eo_flip_used:
                soc_part_m = np.zeros((5, 5)) * 1j
                soc_part_m[3:, :3] = self.sz * self.lattice_params.lamb_m * np.array(
                    [[np.sqrt(3) / 2, -1 / 2, 1j / 2],
                     [-1j / 2 * np.sqrt(3), -1j / 2, -1 / 2]])
                soc_part_m[:3, 3:] = -soc_part_m[3:, :3].T
                reorder_keys = [np.abs(np.array([0, 2, -2, 1, -1]) - key).argmin()
                                for key in self.orbital.l[self.orb_type(self.m_name)]]
                soc_part_m = self._reorder(soc_part_m, [reorder_keys] * 2)
                h_0_m[:5, 5:] = soc_part_m
                h_0_m[5:, :5] = soc_part_m.conj().T
            if self.single_orbital:
                n_m = len(m_orbs)
                for i_m in range(n_m):
                    lat.add_one_sublattice(m_orbs[i_m], [0, 0], np.real(h_0_m[i_m, i_m]))
                for i_m in range(n_m):
                    for j_m in np.arange(i_m + 1, n_m):
                        h_name = self._make_name("h_0_m", 0, m_orbs[i_m], m_orbs[j_m])
                        lat.register_hopping_energies(dict([(h_name, h_0_m[i_m, j_m])]))
                        lat.add_one_hopping([0, 0], m_orbs[i_m], m_orbs[j_m], h_name)
            else:
                lat.add_one_sublattice(self.m_name, [0, 0], h_0_m.T)
            if self.lat4:
                if self.single_orbital:
                    for i_m in range(n_m):
                        lat.add_one_sublattice(m2_orbs[i_m], [0, 0], np.real(h_0_m[i_m, i_m]))
                    for i_m in range(n_m):
                        for j_m in np.arange(i_m + 1, n_m):
                            h_name = self._make_name("h_0_m", 0, m2_orbs[i_m], m2_orbs[j_m])
                            lat.register_hopping_energies(dict([(h_name, h_0_m[i_m, j_m])]))
                            lat.add_one_hopping([0, 0], m2_orbs[i_m], m2_orbs[j_m], h_name)
                else:
                    lat.add_one_sublattice(self.m_name + "2",
                                           [self.lattice_params.a / 2, self.lattice_params.a * np.sqrt(3) / 2],
                                           h_0_m.T)

        if self.lattice_params.h_0_c is not None:
            h_0_c = self._make_onsite(self.lattice_params.h_0_c, self.x_name, self.lattice_params.lamb_c)
            c_orbs = self._c_orbs
            c2_orbs = self._c2_orbs
            if self.soc_eo_flip_used:
                soc_part_c = np.zeros((6, 6)) * 1j
                soc_part_c[:3, 3:] = self.sz * self.lattice_params.lamb_c * np.array(
                    [[0, 0, 1 / 2],
                     [0, 0, -1j / 2],
                     [-1 / 2, 1j / 2, 0]])
                soc_part_c[3:, :3] = -soc_part_c[:3, 3:].T
                reorder_keys1 = [np.abs(np.array([1, -1, 0]) - key).argmin()
                                 for key in self.orbital.l[self.orb_type(self.x_name)][:3]]
                reorder_keys2 = [3 + np.abs(np.array([1, -1, 0]) - key).argmin()
                                 for key in self.orbital.l[self.orb_type(self.x_name)][3:]]
                soc_part_c = self._reorder(soc_part_c, [reorder_keys1 + reorder_keys2] * 2)
                h_0_c[:6, 6:] = soc_part_c
                h_0_c[6:, :6] = soc_part_c.conj().T
            if self.single_orbital:
                n_c = len(c_orbs)
                for i_c in range(n_c):
                    lat.add_one_sublattice(c_orbs[i_c],
                                           [self.lattice_params.a / 2, self.lattice_params.a * np.sqrt(3) / 6],
                                           np.real(h_0_c[i_c, i_c]))
                for i_c in range(n_c):
                    for j_c in np.arange(i_c + 1, n_c):
                        h_name = self._make_name("h_0_c", 0, c_orbs[i_c], c_orbs[j_c])
                        lat.register_hopping_energies(dict([(h_name, h_0_c[i_c, j_c])]))
                        lat.add_one_hopping([0, 0], c_orbs[i_c], c_orbs[j_c], h_name)
            else:
                lat.add_one_sublattice(self.x_name,
                                       [self.lattice_params.a / 2, self.lattice_params.a * np.sqrt(3) / 6],
                                       h_0_c.T)

            if self.lat4:
                if self.single_orbital:
                    for i_c in range(n_c):
                        lat.add_one_sublattice(c2_orbs[i_c],
                                               [self.lattice_params.a / 2, self.lattice_params.a * np.sqrt(3) / 6],
                                               np.real(h_0_c[i_c, i_c]))
                    for i_c in range(n_c):
                        for j_c in np.arange(i_c + 1, n_c):
                            h_name = self._make_name("h_0_c", 0, c2_orbs[i_c], c2_orbs[j_c])
                            lat.register_hopping_energies(dict([(h_name, h_0_c[i_c, j_c])]))
                            lat.add_one_hopping([0, 0], c2_orbs[i_c], c2_orbs[j_c], h_name)
                else:
                    lat.add_one_sublattice(self.x_name + "2",
                                           [0, self.lattice_params.a * 2 * np.sqrt(3) / 3],
                                           h_0_c.T)

        if self.lattice_params.h_1_m is not None:
            cos = [[-1, -1], [0, 0], [-1, 0]] if not self.lat4 else [[0, -1], [0, 0], [0, 0], [1, 0], [-1, 0], [0, 0]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c2, c1, c1, c2, c1, c2] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name + "2", self.x_name, self.x_name, self.x_name + "2", self.x_name, self.x_name + "2"]
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_1_m,
                                    f_n=self.m_name,
                                    t_n=self.x_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_1_m")

        if self.lattice_params.h_2_m is not None:
            cos = [[1, 0], [0, 1], [-1, -1]] if not self.lat4 else [[1, 0], [1, 0], [-1, 0], [0, 1], [-1, -1], [0, 0]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [[m1, m2, m2, m1, m2, m1] for (m1, m2) in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2", self.m_name + "2", self.m_name, self.m_name + "2", self.m_name]
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_2_m,
                                    f_n=self.m_name,
                                    t_n=self.m_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_2_m")

        if self.lattice_params.h_2_c is not None:
            cos = [[1, 0], [0, 1], [-1, -1]] if not self.lat4 else [[1, 0], [1, 0], [0, 0], [-1, 1], [0, -1], [-1, 0]]
            fnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [ci * 3 for ci in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c1, c2, c2, c1, c2, c1] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2", self.x_name + "2", self.x_name, self.x_name + "2", self.x_name]
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_2_c,
                                    f_n=self.x_name,
                                    t_n=self.x_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_2_c")

        if self.lattice_params.h_3_m is not None:
            cos = [[0, 1], [-2, -1], [0, -1]] if not self.lat4 else [[0, 0], [0, 1], [-1, -1], [-1, 0], [1, -1], [1, 0]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [ci * 3 for ci in zip(c2_orbs, c_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name + "2", self.x_name] * 3
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_3_m,
                                    f_n=self.m_name,
                                    t_n=self.x_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_3_m")

        if self.lattice_params.h_4_m is not None:
            _, h_4_ma, h_4_mb = self._make_h_angle(self.lattice_params.h_4_m, self.m_name, self.x_name, np.arctan(np.sqrt(3) / 5))
            cosa = [[-1, -2], [1, 1], [-2, 0]] if not self.lat4 else [[0, -1], [1, -1], [1, 0], [1, 1], [-2, 0], [-1, 0]]
            cosb = [[-2, -2], [1, 0], [-1, 1]] if not self.lat4 else [[-1, -1], [0, -1], [1, 0], [2, 0], [-1, 0], [-1, 1]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnla = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c1, c2, c2, c1, c1, c2] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2", self.x_name + "2", self.x_name, self.x_name, self.x_name + "2"] * 3
                )
            )
            tnlb = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c1, c2, c1, c2, c2, c1] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2", self.x_name, self.x_name + "2", self.x_name + "2", self.x_name] * 3
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=h_4_ma,
                                    f_n=self.m_name,
                                    t_n=self.x_name,
                                    cos=cosa,
                                    fnl=fnl,
                                    tnl=tnla,
                                    h_name="h_4_ma")
            lat = self._add_hopping(lat=lat,
                                    mat=h_4_mb,
                                    f_n=self.m_name,
                                    t_n=self.x_name,
                                    cos=cosb,
                                    fnl=fnl,
                                    tnl=tnlb,
                                    h_name="h_4_mb")

        if self.lattice_params.h_5_m is not None:
            cos = [[1, 2], [-2, -1], [1, -1]] if not self.lat4 else [[0, 1], [0, 1], [-2, -1], [-1, 0], [1, -1], [2, 0]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [[m1, m2, m2, m1, m2, m1] for (m1, m2) in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2", self.m_name + "2", self.m_name, self.m_name + "2", self.m_name]
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_5_m,
                                    f_n=self.m_name,
                                    t_n=self.m_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_5_m")

        if self.lattice_params.h_5_c is not None:
            cos = [[1, 2], [-2, -1], [1, -1]] if not self.lat4 else [[0, 1], [0, 1], [-1, -1], [-2, 0], [2, -1], [1, 0]]
            fnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [ci * 3 for ci in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c1, c2, c2, c1, c2, c1] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2", self.x_name + "2", self.x_name, self.x_name + "2", self.x_name]
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_5_c,
                                    f_n=self.x_name,
                                    t_n=self.x_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_5_c")

        if self.lattice_params.h_6_m is not None:
            cos = [[2, 0], [0, 2], [-2, -2]] if not self.lat4 else [[2, 0], [2, 0], [-1, 1], [-1, 1], [-1, -1], [-1, -1]]
            fnl = (
                (
                    [[m_i] * 3 for m_i in m_orbs]
                    if not self.lat4 else
                    [mi * 3 for mi in zip(m_orbs, m2_orbs)]
                ) if self.single_orbital else (
                    [self.m_name] * 3
                    if not self.lat4 else
                    [self.m_name, self.m_name + "2"] * 3
                )
            )
            tnl = fnl
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_6_m,
                                    f_n=self.m_name,
                                    t_n=self.m_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_6_m")

        if self.lattice_params.h_6_c is not None:
            cos = [[2, 0], [0, 2], [-2, -2]] if not self.lat4 else [[2, 0], [2, 0], [-1, 1], [-1, 1], [-1, -1], [-1, -1]]
            fnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [ci * 3 for ci in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2"] * 3
                )
            )
            tnl = (
                (
                    [[c_i] * 3 for c_i in c_orbs]
                    if not self.lat4 else
                    [[c1, c2, c2, c1, c2, c1] for (c1, c2) in zip(c_orbs, c2_orbs)]
                ) if self.single_orbital else (
                    [self.x_name] * 3
                    if not self.lat4 else
                    [self.x_name, self.x_name + "2"] * 3
                )
            )
            lat = self._add_hopping(lat=lat,
                                    mat=self.lattice_params.h_6_c,
                                    f_n=self.x_name,
                                    t_n=self.x_name,
                                    cos=cos,
                                    fnl=fnl,
                                    tnl=tnl,
                                    h_name="h_6_c")
        return lat
