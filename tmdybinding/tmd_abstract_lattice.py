"""Abstract classes for a tmd-lattice, including for the storage of lattice variables and the orbitals."""
import pybinding as pb
import numpy as np
import re
from typing import Optional, List, Tuple, Dict
from .parameters import ParametersList
from abc import ABC, abstractmethod


class VariableStorage:
    """
    Class for saving the matrices and variables in the AbstractLattice-class.
    There are also checks for a consistent model, where the sizes of the matrices, the included hoppings and the onsite
    matrices are checked.
    """

    def __init__(self, params_dict: Optional[dict] = None):
        """Make a VariableStorage object with the given parameters. If no parameters are given, the parameters are None.

        Parameters:
            params_dict (Optional[dict]): The parameters for the lattice model.
                The dictionary entries must be in the format `"name": parameter`.
                The possible parameters are `"h_0_m"`, `"h_0_c"`, `"h_1_m"`, `"h_2_m"`, `"h_2_c"`, `"h_3_m"`, `"h_4_m"`,
                `"h_5_m"`, `"h_5_c"`, `"h_6_m"`, `"h_6_c"`, `"a"`, `"lamb_m"` and `"lamb_c"`.
        """
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
        if params_dict is not None:
            self.set_params(params_dict)

    @staticmethod
    def _params_check(params_dict: dict, param_name: str) -> Optional[np.ndarray]:
        if param_name in params_dict.keys():
            if params_dict[param_name] is not None:
                return np.array(params_dict[param_name])
        return None

    def set_params(self, params_dict: dict, keep: bool = True):
        """Set the parameters for the lattice model.

        Parameters:
            params_dict (dict): The parameters for the lattice model.
                The dictionary entries must be in the format `"name": parameter`.
                The possible parameters are `"h_0_m"`, `"h_0_c"`, `"h_1_m"`, `"h_2_m"`, `"h_2_c"`, `"h_3_m"`, `"h_4_m"`,
                `"h_5_m"`, `"h_5_c"`, `"h_6_m"`, `"h_6_c"`, `"a"`, `"lamb_m"` and `"lamb_c"`.
            keep (bool): If True, the parameters that are not included in the `param_dict`are kept as they are.
                If False, the parameters that are not given are set to None.
        """
        params = self.to_dict() if keep else {}
        for name in params_dict.keys():
            assert str(name) in self.__keys, f"key {name} not in expected hoppings, possible wrong name"
            params[name] = params_dict[name]
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

    def _make_bools(self, params_dict: dict):
        return [self._attr_check(name, params_dict) for name in self.__keys]

    def _component_check(self, params_dict: dict):
        (m_bool, c_bool,
         h_1_m_bool,
         h_2_m_bool, h_2_c_bool,
         h_3_m_bool,
         h_4_m_bool,
         h_5_m_bool, h_5_c_bool,
         h_6_m_bool, h_6_c_bool) = self._make_bools(params_dict)[:-3]
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
        assert h_1_m_bool or h_2_c_bool or h_2_m_bool or h_3_m_bool or h_4_m_bool or h_5_c_bool or h_5_m_bool or \
               h_6_c_bool or h_6_m_bool, "no hoppings specified in the model"

    def _shape_check(self, params_dict: dict):
        (m_bool, c_bool,
         h_1_m_bool,
         h_2_m_bool, h_2_c_bool,
         h_3_m_bool,
         h_4_m_bool,
         h_5_m_bool, h_5_c_bool,
         h_6_m_bool, h_6_c_bool) = self._make_bools(params_dict)[:-3]
        if m_bool:
            m_shape = self._check_shape(params_dict["h_0_m"])
            if h_1_m_bool:
                assert np.shape(params_dict["h_1_m"])[1] == m_shape, "shape 1st hopping from M not correct"
            if h_2_m_bool:
                assert np.shape(params_dict["h_2_m"])[0] == m_shape, "shape 2nd hopping to M not correct"
                assert np.shape(params_dict["h_2_m"])[1] == m_shape, "shape 2nd hopping from M not correct"
            if h_3_m_bool:
                assert np.shape(params_dict["h_3_m"])[1] == m_shape, "shape 3rd hopping from M not correct"
            if h_4_m_bool:
                assert np.shape(params_dict["h_4_m"])[1] == m_shape, "shape 4rd hopping from M not correct"
            if h_5_m_bool:
                assert np.shape(params_dict["h_5_m"])[0] == m_shape, "shape 5th hopping to M not correct"
                assert np.shape(params_dict["h_5_m"])[1] == m_shape, "shape 5th hopping from M not correct"
            if h_6_m_bool:
                assert np.shape(params_dict["h_6_m"])[0] == m_shape, "shape 6th hopping to M not correct"
                assert np.shape(params_dict["h_6_m"])[1] == m_shape, "shape 6th hopping from M not correct"
        if c_bool:
            c_shape = self._check_shape(params_dict["h_0_c"])
            if h_1_m_bool:
                assert np.shape(params_dict["h_1_m"])[0] == c_shape, "shape 1st hopping to X not correct"
            if h_2_c_bool:
                assert np.shape(params_dict["h_2_c"])[0] == c_shape, "shape 2nd hopping to X not correct"
                assert np.shape(params_dict["h_2_c"])[1] == c_shape
            if h_3_m_bool:
                assert np.shape(params_dict["h_3_m"])[0] == c_shape, "shape 3rd hopping to X not correct"
            if h_4_m_bool:
                assert np.shape(params_dict["h_4_m"])[0] == c_shape, "shape 4rd hopping to X not correct"
            if h_5_c_bool:
                assert np.shape(params_dict["h_5_c"])[0] == c_shape, "shape 5th hopping to X not correct"
                assert np.shape(params_dict["h_5_c"])[1] == c_shape, "shape 5th hopping from X not correct"
            if h_6_c_bool:
                assert np.shape(params_dict["h_6_c"])[0] == c_shape, "shape 6th hopping to X not correct"
                assert np.shape(params_dict["h_6_c"])[1] == c_shape, "shape 6th hopping from X not correct"

    @staticmethod
    def _attr_check(name, params_dict: dict):
        if name in params_dict.keys():
            return params_dict[name] is not None
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

    def to_dict(self) -> dict:
        """Return the parameters as a dictionary in the format `"name": parameter`."""
        dict_list = []
        for name in self.__keys:
            value = getattr(self, name)
            if value is not None:
                dict_list.append((name, value))
        return dict(dict_list)

    @property
    def h_0_m(self) -> np.ndarray:
        """The onsite energy for the metal atom."""
        return self.__h_0_m

    @h_0_m.setter
    def h_0_m(self, value: np.ndarray):
        self.set_params({"h_0_m": value})

    @property
    def h_0_c(self) -> np.ndarray:
        """The onsite energy for the chalcogen atom."""
        return self.__h_0_c

    @h_0_c.setter
    def h_0_c(self, value: np.ndarray):
        self.set_params({"h_0_c": value})

    @property
    def h_1_m(self) -> np.ndarray:
        """The first nearest neighbour hopping from the metal atom."""
        return self.__h_1_m

    @h_1_m.setter
    def h_1_m(self, value: np.ndarray):
        self.set_params({"h_1_m": value})

    @property
    def h_2_m(self) -> np.ndarray:
        """The second-nearest neighbour hopping from the metal atom."""
        return self.__h_2_m

    @h_2_m.setter
    def h_2_m(self, value: np.ndarray):
        self.set_params({"h_2_m": value})

    @property
    def h_2_c(self) -> np.ndarray:
        """The second-nearest neighbour hopping from the chalcogen atom."""
        return self.__h_2_c

    @h_2_c.setter
    def h_2_c(self, value: np.ndarray):
        self.set_params({"h_2_c": value})

    @property
    def h_3_m(self) -> np.ndarray:
        """The third-nearest neighbour hopping from the metal atom."""
        return self.__h_3_m

    @h_3_m.setter
    def h_3_m(self, value: np.ndarray):
        self.set_params({"h_3_m": value})

    @property
    def h_4_m(self) -> np.ndarray:
        """The fourth-nearest neighbour hopping from the metal atom."""
        return self.__h_4_m

    @h_4_m.setter
    def h_4_m(self, value: np.ndarray):
        self.set_params({"h_4_m": value})

    @property
    def h_5_m(self) -> np.ndarray:
        """The fifth-nearest neighbour hopping from the metal atom."""
        return self.__h_5_m

    @h_5_m.setter
    def h_5_m(self, value: np.ndarray):
        self.set_params({"h_5_m": value})

    @property
    def h_5_c(self) -> np.ndarray:
        """The fifth-nearest neighbour hopping from the chalcogen atom."""
        return self.__h_5_c

    @h_5_c.setter
    def h_5_c(self, value: np.ndarray):
        self.set_params({"h_5_c": value})

    @property
    def h_6_m(self) -> np.ndarray:
        """The sixth-nearest neighbour hopping from the metal atom."""
        return self.__h_6_m

    @h_6_m.setter
    def h_6_m(self, value: np.ndarray):
        self.set_params({"h_6_m": value})

    @property
    def h_6_c(self) -> np.ndarray:
        """The sixth-nearest neighbour hopping from the chalcogen atom."""
        return self.__h_6_c

    @h_6_c.setter
    def h_6_c(self, value: np.ndarray):
        self.set_params({"h_6_c": value})

    @property
    def a(self) -> float:
        """The lattice constant."""
        return self.__a

    @a.setter
    def a(self, value: float):
        self.set_params({"a": value})

    @property
    def lamb_m(self) -> float:
        """The spin-orbit coupling for the metal atom."""
        return self.__lamb_m

    @lamb_m.setter
    def lamb_m(self, value: float):
        self.set_params({"lamb_m": value})

    @property
    def lamb_c(self) -> float:
        """The spin-orbit coupling for the chalcogen atom."""
        return self.__lamb_c

    @lamb_c.setter
    def lamb_c(self, value: float):
        self.set_params({"lamb_c": value})


class LatticeOrbitals:
    """ Object for saving the l and s values of a lattice model."""

    def __init__(self, l_number: Dict[str, List[int]], orbs: Dict[str, List[str]],
                 group: Optional[Dict[str, List[int]]] = None, clockwise: bool = False):
        """Make a LatticeOrbitals object with the given combinations of l and orbs.
        The input for `l`, `orbs` and `group` must be in the format `"atom_name": variables`.


        Parameters:
            l_number (Dict[str, List[int]]): The l-values for the orbitals.
                The dictionary entries must be in the format `"atom_name": l_number`.
                The l-values must be integers.
                They must come in pairs of negative and positive values, except for zero.
                Example: `{"S": [0, 1, -1, 0, 1, -1]}
            orbs (Dict[str, List[str]]): The orbital names for the orbitals.
                The dictionary entries must be in the format `"atom_name": orbs`.
                The orbital names must be strings.
                Example: `{"S": ["pze", "pxe", pye", "pzo", "pxo", "pyo"]}`
            group (Optional[Dict[str, List[int]]): The group of orbitals with the same l-value.
                If multiple orbitals with the same `l` are present on one atom, you can group the similar orbitals.
                If all the `l`-values are unique, the group is not needed.
                The dictionary entries must be in the format `"atom_name": group`.
                Example: `{"S": [0, 1, 1, 2, 3, 3]}`
            clockwise (bool): If True, the rotation matrices are clockwise.
        """
        self.__l_number = None
        self.__orbs = None
        self.__names = None
        self.__group = None
        self.clockwise = clockwise
        self.set_params(l_number=l_number, orbs=orbs, group=group)

    def set_params(self, l_number: Optional[Dict[str, List[int]]] = None, orbs: Optional[Dict[str, List[str]]] = None,
                   group: Optional[Dict[str, List[int]]] = None):
        """Set the parameters for the lattice model.

        Parameters:
            l_number (Optional[Dict[str, List[int]]]): The l-values for the orbitals.
                The dictionary entries must be in the format `"atom_name": l_number`.
                The l-values must be integers.
                They must come in pairs of negative and positive values, except for zero.
                If None, the previous value is kept.
                Example: `{"S": [0, 1, -1, 0, 1, -1]}
            orbs (Optional[Dict[str, Dict[str]]]): The orbital names for the orbitals.
                The dictionary entries must be in the format `"atom_name": orbs`.
                The orbital names must be strings.
                If None, the previous value is kept.
                Example: `{"S": ["pze", "pxe", pye", "pzo", "pxo", "pyo"]}`
            group (Optional[Dict[str, List[int]]): The group of orbitals with the same l-value.
                If multiple orbitals with the same `l` are present on one atom, you can group the similar orbitals.
                If all the `l`-values are unique, the group is not needed.
                If None, the previous value is kept (if there).
                The dictionary entries must be in the format `"atom_name": group`.
                Example: `{"S": [0, 1, 1, 2, 3, 3]}`
        """
        l_pass = l_number or self.l_number
        orbs_pass = orbs or self.orbs
        group_pass = group or self.group
        self._check_names(l_pass, orbs_pass, group_pass)
        self._check_shape(l_pass, orbs_pass, group_pass)

    @staticmethod
    def _make_bools(l_number: Dict[str, List[int]], orbs: Dict[str, List[str]],
                    group: Optional[Dict[str, List[int]]]) -> Tuple[bool, bool, bool]:
        """Make the boolean values for the input.

        Parameters:
            l_number (Dict[str, List[int]]): The l-values for the orbitals.
            orbs (Dict[str, List[str]]): The orbital names for the orbitals.
            group (Optional[Dict[str, List[int]]): The group of orbitals with the same l-value.
        """
        l_bool = l_number is not None
        orbs_bool = orbs is not None
        group_bool = group is not None
        assert l_bool or orbs_bool, "not enough input given, at least l or orbs needed"
        return l_bool, orbs_bool, group_bool

    def _check_names(self, l_number: Dict[str, List[int]], orbs: Dict[str, List[str]],
                     group: Optional[Dict[str, List[int]]]):
        """Check if the names of the orbitals are the same."""
        (l_bool, orbs_bool, group_bool) = self._make_bools(l_number, orbs, group)
        names = None
        if l_bool:
            names_l = [str(name) for name in l_number.keys()]
            names = names_l
        if orbs_bool:
            names_orbs = [str(name) for name in orbs.keys()]
            if names is None:
                names = names_orbs
            else:
                assert names == names_orbs, f"the names of l_number and orbs are not the same, {names} != {names_orbs}"
        if group_bool:
            names_group = [str(name) for name in group.keys()]
            assert names == names_group, \
                f"the names of l_number/orbs  and group are not the same, {names} != {names_group}"
        self.__names = names

    def _check_shape(self, l_number: Dict[str, List[int]], orbs: Dict[str, List[str]],
                     group: Optional[Dict[str, List[int]]]):
        """Check if the shape of the orbitals is correct."""
        (l_bool, orbs_bool, group_bool) = self._make_bools(l_number, orbs, group)
        shape = None
        if l_bool:
            shape_l = [np.shape(l_number[name]) for name in self.names]
            shape = shape_l
        if orbs_bool:
            shape_orbs = [np.shape(orbs[name]) for name in self.names]
            if shape is None:
                shape = shape_orbs
            else:
                assert shape == shape_orbs, f"the shape of l_number and orbs are not the same, {shape} !m {shape_orbs}"
        if group_bool:
            shape_group = [np.shape(group[name]) for name in self.names]
            assert shape == shape_group, "the shape of group and l_number and/or orbs are not the same"

        # see if definition of group and l are correct. If no group is specified, but there is an l, the l is checked
        l_local = l_number if l_bool else dict([(name, np.zeros(shape[j])) for j, name in enumerate(self.names)])
        if not group_bool:
            for name in self.names:
                l_num = np.array(l_local[name])
                for lm in set(np.abs(l_num)):
                    l_b = np.abs(l_num) == lm
                    assert np.sum(l_b) < 3, \
                        f"the representation can only be given for max two parts, {np.sum(l_b)} given"
                    if np.sum(l_b) == 2:
                        assert l_num[l_b][0] == -l_num[l_b][1], \
                            f"can't have same l-number if group is not defined, for '{lm}'"
                    else:
                        assert lm == 0, \
                            f"can't have a sole l-number other than zero, '{lm}' given"
                    assert np.sum(l_b) == 1, \
                        f"can't have a sole l-number, only one '{lm}' given"
        group_local = group if group_bool else dict([(name, np.abs(l_local[name])) for name in self.names])

        # check l and group
        for name in self.names:
            group_l = np.array(group_local[name])
            l_num = np.array(l_local[name])
            for group_i in set(group_l):
                group_b = group_i == group_l
                assert np.sum(group_b) < 3, "the group for representation can only be given for max two parts"
                assert (np.sum(group_b) == 2 and l_num[group_b][0] == -l_num[group_b][1]) or np.sum(group_b) == 1, \
                    "can't have same l-number in same group"
        self.__l_number = l_local
        self.__orbs = orbs
        self.__group = group_local

    @staticmethod
    def rot_mat(phi: float = 0) -> np.ndarray:
        """Make a rotation matrix for the given angle.

        Parameters:
            phi (float): The angle of the rotation matrix.
        """
        return np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

    @property
    def l_number(self) -> Dict[str, List[int]]:
        """The l-values for the orbitals."""
        return self.__l_number

    @property
    def names(self) -> List[str]:
        """The names of the atoms."""
        return self.__names

    @property
    def orbs(self) -> Dict[str, List[str]]:
        """The orbital names for the orbitals."""
        return self.__orbs

    @property
    def group(self) -> Dict[str, List[int]]:
        """The group of orbitals with the same l-value."""
        return self.__group

    def _make_matrix(self, matrix_func, single) -> Dict[str, np.ndarray]:
        out_dict: Dict[str, np.ndarray] = {}
        for name in self.names:
            group_l = np.array(self.group[name])
            l_num = np.array(self.l_number[name])
            value = np.zeros((len(l_num), len(l_num)))
            for group_i in set(group_l):
                group_b = group_i == group_l
                lm = np.abs(l_num[group_b][0])
                sign = -1 if l_num[group_b][0] < 0 else 1
                if np.sum(group_b) == 2:
                    matrix = matrix_func(lm, sign)
                    for ik, k in enumerate(np.arange(len(l_num))[group_b]):
                        for ikk, kk in enumerate(np.arange(len(l_num))[group_b]):
                            value[k, kk] = matrix[ik, ikk]
                else:
                    value[np.arange(len(l_num))[group_b], np.arange(len(l_num))[group_b]] = single(lm, sign)
            out_dict[name] = value
        return out_dict

    @property
    def ur(self) -> Dict[str, np.ndarray]:
        """Rotation matrix for the orbitals over an angle of 2pi/3.
        Returns a dict of rotation matrices for the orbitals in the format `"atom_name": matrix`.
        """

        def matrix_func(lm, sign):
            ur_t = self.rot_mat(np.pi * 2 / 3 * lm)
            if sign == -1:
                ur_t = np.transpose(ur_t)
            if self.clockwise:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 1 + (lm * 0 + sign * 0)

        return self._make_matrix(matrix_func, single)

    @property
    def sr(self) -> Dict[str, np.ndarray]:
        """Mirror on yz-plane.
        Returns a dict of mirror matrices for the orbitals in the format `"atom_name": matrix`.
        """

        def matrix_func(lm, sign):
            return sign * np.diag([-1, 1]) * (1 if np.abs(lm) == 1 else -1)

        def single(lm, sign):
            return 1 + (lm * 0 + sign * 0)

        return self._make_matrix(matrix_func, single)

    def ur_angle(self, angle) -> Dict[str, np.ndarray]:
        """Rotation matrix for the orbitals over an angle of `angle`.
        Returns a dict of rotation matrices for the orbitals over an angle of `angle` in the
        format `"atom_name": matrix`.

        Parameters:
            angle (float): The angle of the rotation matrix.
        """

        def matrix_func(lm, sign):
            ur_t = self.rot_mat(angle * lm)
            if sign == -1:
                ur_t = np.transpose(ur_t)
            if self.clockwise:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 1 + (lm * 0 + sign * 0)

        return self._make_matrix(matrix_func, single)

    @property
    def s_h(self):
        """Spin factor for the orbitals.
        Returns a dict of spin matrices for the orbitals in the format `"atom_name": matrix`.
        """

        def matrix_func(lm, sign):
            ur_t = lm / 2 * np.array([[0, -1], [1, 0]])
            if sign == -1:
                ur_t = np.transpose(ur_t)
            return ur_t

        def single(lm, sign):
            return 0 + (lm * 0 + sign * 0)

        return self._make_matrix(matrix_func, single)


class AbstractLattice(ABC):
    """Abstract class for a tmd-lattice"""

    def __init__(self, orbital: LatticeOrbitals, params: ParametersList, soc: bool = False,
                 soc_eo_flip: bool = False, lat4: bool = False, single_orbital: bool = False, soc_sz_part: bool = True,
                 soc_polarized: bool = False, soc_sz: float = 1.,
                 n_v: int = 0, n_b: int = 0, lattice_name: str = "Abstract Lattice"):
        """Make an AbstractLattice object with the given parameters.

        Parameters:
            orbital (LatticeOrbitals): The orbitals for the lattice model.
            params (ParametersList): The parameters for the lattice model.
            soc (bool): If True, the spin-orbit coupling is considere.
            soc_eo_flip (bool): If True, the spin-flip term is included.
            lat4 (bool): If True, the lattice is a 4-atom lattice for use with the armchair direction.
            single_orbital (bool): If True, the lattice is a combination of suborbitals (slower as no matrices are used)
            soc_sz_part (bool): If False, the hamiltonian will not include the Sz part of the spin-orbit coupling.
            soc_polarized (bool): If True, the spin-orbit coupling is polarized (only spin-up part, given by `soc_sz`).
            soc_sz (float): The Sz part of the spin-orbit coupling.
            n_v (int): The index of the highest valence band (in the unitcell).
            n_b (int): The number of bands (in the unitcell).
            lattice_name (str): The name of the lattice. The name of the matrial is obtained from `params`.
        """
        self.lattice_params = VariableStorage()
        self.orbital: LatticeOrbitals = orbital
        self.__n_valence_band: int = n_v
        self.__n_bands: int = n_b
        self.__name: str = "MoS2"
        self.__x_name: str = "S"
        self.__x_type: str = "X"
        self.__m_name: str = "Mo"
        self.__m_type: str = "M"
        self.__params: ParametersList = ParametersList()
        self.__soc_eo_flip: bool = False
        self.__lattice_name: str = lattice_name
        self._lat4: bool = lat4
        self.single_orbital: bool = single_orbital
        self.soc: bool = soc
        self.soc_polarized: bool = soc_polarized
        self.soc_sz_part: bool = soc_sz_part
        self.sz: float = soc_sz
        self.params = params
        self.soc_eo_flip = soc_eo_flip

    @property
    def soc_eo_flip(self) -> bool:
        """If True, the spin-flip term is included.
        The structure is checked for the right shape and properties that it supports a spin-flip term.
        """
        return self.__soc_eo_flip

    @soc_eo_flip.setter
    def soc_eo_flip(self, pol_bool: bool):
        if pol_bool:
            if self.lattice_params.h_0_m is not None:
                assert len(self.orbital.l_number[self.orb_type(self.m_name)]) == 5, \
                    "the metal doesn't have the right shape for spin-flip term"
                assert sorted(self.orbital.l_number[self.orb_type(self.m_name)]) == sorted([0, 2, -2, 1, -1]), \
                    "the metal l is wrong for spin-flip"
            if self.lattice_params.h_0_c is not None:
                assert len(self.orbital.l_number[self.orb_type(self.x_name)]) == 6, \
                    "the chal. doesn't have the right shape for spin-flip term"
                assert sorted(self.orbital.l_number[self.orb_type(self.x_name)]) == sorted([1, -1, 0, 1, -1, 0]), \
                    "the chal. l is wrong for spin-flip"
        self.__soc_eo_flip = pol_bool

    @property
    def params(self) -> ParametersList:
        """The parameters for the lattice model."""
        return self.__params

    @params.setter
    def params(self, params: ParametersList):
        self.__params = params
        self.name = self.params["material"]

    @property
    def lattice_name(self) -> str:
        """The name of the lattice. This is not the name of the material."""
        return self.__lattice_name

    @lattice_name.setter
    def lattice_name(self, lattice_name: str):
        self.__lattice_name = lattice_name

    @property
    def name(self) -> str:
        """The name of the material."""
        return self.__name

    @property
    def m_name(self) -> str:
        """The name of the metal atom."""
        return self.__m_name

    def orb_type(self, z_name: str) -> Optional[str]:
        """Return the orbital type of the atom with the given name. Usefull for the case when `lat4` is True.
        Returns `"M"` for the metal atom and `"X"` for the chalcogen atom.
        """
        if z_name == self.m_name or (self.lat4 and z_name == self.m_name + "2"):
            return self.__m_type
        elif z_name == self.x_name or (self.lat4 and z_name == self.x_name + "2"):
            return self.__x_type
        else:
            return None

    def z_name(self, orb_type: str) -> Optional[str]:
        """Return the name of the atom with the given orbital type, for example, obtain `"S"` from `"X"` in the case
        the material is `"MoS2"`.
        """
        if orb_type == self.__m_type:
            return self.m_name
        elif orb_type == self.__x_type:
            return self.x_name
        else:
            return None

    @property
    def x_name(self) -> str:
        """The name of the chalcogen atom."""
        return self.__x_name

    @name.setter
    def name(self, name: str):
        self.__name = name
        self.__m_name, self.__x_name = re.findall("[A-Z][a-z]*", self.name)
        self._generate_matrices()

    @property
    def lat4(self) -> bool:
        """If True, the lattice is a 4-atom lattice for use with the armchair direction."""
        return self._lat4

    @lat4.setter
    def lat4(self, lat4: bool):
        self._lat4 = lat4
        self._generate_matrices()

    @property
    def a1(self) -> np.ndarray:
        """The first lattice vector. Corrected if `lat4` is True."""
        return np.array([1, 0]) * self.lattice_params.a

    @property
    def a2(self) -> np.ndarray:
        """The second lattice vector. Corrected if `lat4` is True."""
        return (np.array([0, np.sqrt(3)]) if self.lat4 else np.array([-1 / 2, np.sqrt(3) / 2])) * self.lattice_params.a

    @property
    def soc_doubled_ham(self) -> bool:
        """If True, the hamiltonian is doubled for the spin-orbit coupling."""
        return not self.soc_polarized if self.soc else False

    @property
    def soc_eo_flip_used(self) -> bool:
        """If True, the spin-flip term is used."""
        return self.soc and self.soc_eo_flip and not self.soc_polarized

    @property
    def n_valence_band(self) -> int:
        """The index of the highest valence band. Corrected for the SOC and `lat4`."""
        return (self.__n_valence_band + 1) * 2 - 1 if self.soc_doubled_ham else self.__n_valence_band

    @property
    def n_bands(self) -> int:
        """The total number of bands. Corrected for the SOC and `lat4`."""
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
        return str(out[0]), int(out[1]), str(out[2]), str(out[3])

    @abstractmethod
    def _generate_matrices(self):
        """Generate the matrices for the lattice model.
        This function should be implemented in the child class.
        """
        pass

    @staticmethod
    def _ham(h: np.ndarray, ur_l: np.ndarray, ur_r: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return h, ur_l.dot(h.dot(ur_r.T)), ur_l.T.dot(h.dot(ur_r))

    @staticmethod
    def _reorder(matrix: np.ndarray, keys: Tuple[List[int], List[int]]) -> np.ndarray:
        return np.array([[matrix[xi, yi] for yi in keys[1]] for xi in keys[0]])

    def block_diag(self, *matrices) -> np.ndarray:
        """Make a block diagonal matrix from the given matrices."""
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
            return np.concatenate((np.concatenate((matrix0, np.zeros(z_1)), axis=nd - 1),
                                   np.concatenate((np.zeros(z_2), matrix1), axis=nd - 1)), axis=nd - 2)

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
            hn = [hn[int(ih / 2)] for ih in range(6)]
            n_n_n = [0, 0, 1, 1, 2, 2]
        # n_n = 6 if self.lat4 else 3
        if self.single_orbital:
            for f_i, nf_i in enumerate(fnl):
                for t_j, nt_j in enumerate(tnl):
                    h_names = [self._make_name(h_name, n_i, nfi, ntj) for n_i, nfi, ntj in zip(n_n_n, nf_i, nt_j)]
                    lat.register_hopping_energies(dict([(h_n_i, h.conj().T[f_i, t_j]) for h_n_i, h in zip(h_names, hn)]))
                    lat.add_hoppings(*[(co, nfi, ntj, h_n_i) for co, h_n_i, nfi, ntj in zip(cos, h_names, nf_i, nt_j)])
        else:
            h_names = [self._make_name(h_name, n_i, nfi, ntj) for n_i, nfi, ntj in zip(n_n_n, fnl, tnl)]
            lat.register_hopping_energies(dict([(h_n_i, h.conj().T) for h_n_i, h in zip(h_names, hn)]))
            lat.add_hoppings(*[(co, nfi, ntj, h_n_i) for co, h_n_i, nfi, ntj in zip(cos, h_names, fnl, tnl)])
        return lat

    def lattice(self) -> pb.Lattice:
        """Make the lattice model. It returns a `pybinding.Lattice` object."""
        lat = pb.Lattice(a1=self.a1, a2=self.a2)
        m_orbs, m2_orbs = 0, 0
        c_orbs, c2_orbs = 0, 0
        if self.lattice_params.h_0_m is not None:
            m_orbs = self._m_orbs
            m2_orbs = self._m2_orbs
            h_0_m = self._make_onsite(self.lattice_params.h_0_m, self.m_name, self.lattice_params.lamb_m)
            n_m, n_c = 0, 0
            if self.soc_eo_flip_used:
                soc_part_m = np.zeros((5, 5)) * 1j
                soc_part_m[3:, :3] = self.sz * self.lattice_params.lamb_m * np.array(
                    [[np.sqrt(3) / 2, -1 / 2, 1j / 2],
                     [-1j / 2 * np.sqrt(3), -1j / 2, -1 / 2]])
                soc_part_m[:3, 3:] = -soc_part_m[3:, :3].T
                reorder_keys = [np.abs(np.array([0, 2, -2, 1, -1]) - key).argmin()
                                for key in self.orbital.l_number[self.orb_type(self.m_name)]]
                soc_part_m = self._reorder(soc_part_m, (reorder_keys, reorder_keys))
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
                lat.add_one_sublattice(self.m_name, [0, 0], h_0_m)
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
                                           h_0_m)

        if self.lattice_params.h_0_c is not None:
            h_0_c = self._make_onsite(self.lattice_params.h_0_c, self.x_name, self.lattice_params.lamb_c)
            c_orbs = self._c_orbs
            c2_orbs = self._c2_orbs
            n_c = len(c_orbs)
            if self.soc_eo_flip_used:
                soc_part_c = np.zeros((6, 6)) * 1j
                soc_part_c[:3, 3:] = self.sz * self.lattice_params.lamb_c * np.array(
                    [[0, 0, 1 / 2],
                     [0, 0, -1j / 2],
                     [-1 / 2, 1j / 2, 0]])
                soc_part_c[3:, :3] = -soc_part_c[:3, 3:].T
                reorder_keys1 = [np.abs(np.array([1, -1, 0]) - key).argmin()
                                 for key in self.orbital.l_number[self.orb_type(self.x_name)][:3]]
                reorder_keys2 = [3 + np.abs(np.array([1, -1, 0]) - key).argmin()
                                 for key in self.orbital.l_number[self.orb_type(self.x_name)][3:]]
                soc_part_c = self._reorder(soc_part_c, (reorder_keys1 + reorder_keys2, reorder_keys1 + reorder_keys2))
                h_0_c[:6, 6:] = soc_part_c
                h_0_c[6:, :6] = soc_part_c.conj().T
            if self.single_orbital:
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
                                       h_0_c)

            if self.lat4:
                if self.single_orbital:
                    for i_c in range(n_c):
                        lat.add_one_sublattice(c2_orbs[i_c],
                                               [0, self.lattice_params.a * 2 * np.sqrt(3) / 6],
                                               np.real(h_0_c[i_c, i_c]))
                    for i_c in range(n_c):
                        for j_c in np.arange(i_c + 1, n_c):
                            h_name = self._make_name("h_0_c", 0, c2_orbs[i_c], c2_orbs[j_c])
                            lat.register_hopping_energies(dict([(h_name, h_0_c[i_c, j_c])]))
                            lat.add_one_hopping([0, 0], c2_orbs[i_c], c2_orbs[j_c], h_name)
                else:
                    lat.add_one_sublattice(self.x_name + "2",
                                           [0, self.lattice_params.a * 2 * np.sqrt(3) / 3],
                                           h_0_c)

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
            _, h_4_ma, h_4_mb = self._make_h_angle(
                self.lattice_params.h_4_m, self.m_name, self.x_name, np.arctan(np.sqrt(3) / 5)
            )
            cosa = [[-1, -2], [1, 1], [-2, 0]] \
                if not self.lat4 else [[0, -1], [1, -1], [1, 0], [1, 1], [-2, 0], [-1, 0]]
            cosb = [[-2, -2], [1, 0], [-1, 1]] \
                if not self.lat4 else [[-1, -1], [0, -1], [1, 0], [2, 0], [-1, 0], [-1, 1]]
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
            cos = [[2, 0], [0, 2], [-2, -2]] \
                if not self.lat4 else [[2, 0], [2, 0], [-1, 1], [-1, 1], [-1, -1], [-1, -1]]
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
            cos = [[2, 0], [0, 2], [-2, -2]] \
                if not self.lat4 else [[2, 0], [2, 0], [-1, 1], [-1, 1], [-1, -1], [-1, -1]]
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
