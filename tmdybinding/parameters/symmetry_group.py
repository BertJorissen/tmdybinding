"""The parameters used in the Symmetry-Group models for TMDs."""
from dataclasses import dataclass
from typing import Optional, Dict, Union, List
import warnings


@dataclass
class Parameter:
    """Class to store one separate parameter"""
    name: str = ""


@dataclass
class FloatParameter(Parameter):
    """Class to store one separate float parameter"""
    param: Optional[float] = None


@dataclass
class StringParameter(Parameter):
    """Class to store one separate string parameter"""
    param: Optional[str] = None


class ParametersList:
    """Class to save the parameters"""

    def __init__(self, input_dict: Optional[Dict[str, Union[float, str]]] = None):
        """Initialize the parameters for the TMD lattice.

        Parameters:
            input_dict (Optional[Dict[str, Union[float, str]]]): The parameters for the TMD lattice.
                The keys are the names of the parameters and the values are the values of the parameters.
                The keys are:
                `a`,            `lamb_m`,       `lamb_x`,       `material`,
                `eps_0_x_e`, `eps_1_x_e`, `eps_0_m_e`, `eps_1_m_e`,
                `u_1_0_m_e`, `u_1_1_m_e`, `u_1_2_m_e`, `u_1_3_m_e`, `u_1_4_m_e`,
                `u_2_0_m_e`, `u_2_1_m_e`, `u_2_2_m_e`, `u_2_3_m_e`, `u_2_4_m_e`, `u_2_5_m_e`,
                `u_2_0_x_e`, `u_2_1_x_e`, `u_2_2_x_e`, `u_2_3_x_e`, `u_2_4_x_e`, `u_2_5_x_e`,
                `u_3_0_m_e`, `u_3_1_m_e`, `u_3_2_m_e`, `u_3_3_m_e`, `u_3_4_m_e`,
                `u_4_0_m_e`, `u_4_1_m_e`, `u_4_2_m_e`, `u_4_3_m_e`, `u_4_4_m_e`,
                `u_5_0_m_e`, `u_5_1_m_e`, `u_5_3_m_e`, `u_5_5_m_e`, `u_5_6_m_e`,
                `u_5_0_x_e`, `u_5_2_x_e`, `u_5_3_x_e`, `u_5_5_x_e`, `u_5_6_x_e`,
                `u_6_0_m_e`, `u_6_1_m_e`, `u_6_2_m_e`, `u_6_3_m_e`, `u_6_4_m_e`, `u_6_5_m_e`,
                `u_6_0_x_e`, `u_6_1_x_e`, `u_6_2_x_e`, `u_6_3_x_e`, `u_6_4_x_e`, `u_6_5_x_e`,
                `eps_0_x_o`, `eps_1_x_o`, `eps_0_m_o`,
                `u_1_0_m_o`, `u_1_1_m_o`, `u_1_2_m_o`,
                `u_2_0_m_o`, `u_2_1_m_o`, `u_2_2_m_o`,
                `u_2_0_x_o`, `u_2_1_x_o`, `u_2_2_x_o`, `u_2_3_x_o`, `u_2_4_x_o`, `u_2_5_x_o`,
                `u_3_0_m_o`, `u_3_1_m_o`, `u_3_2_m_o`,
                `u_4_0_m_o`, `u_4_1_m_o`, `u_4_2_m_o`,
                `u_5_0_m_o`, `u_5_2_m_o`,
                `u_5_0_x_o`, `u_5_2_x_o`, `u_5_3_x_o`, `u_5_5_x_o`, `u_5_6_x_o`,
                `u_6_0_m_o`, `u_6_1_m_o`, `u_6_2_m_o`,
                `u_6_0_x_o`, `u_6_1_x_o`, `u_6_2_x_o`, `u_6_3_x_o`, `u_6_4_x_o` and `u_6_5_x_o`
        """
        energy_params = [
            *[f"eps_{i}_x_{r}" for i in range(2) for r in ("e", "o")],
            *[f"eps_{i}_m_{r}" for i in "0" for r in ("e", "o")],
            *[f"eps_{i}_m_{r}" for i in "1" for r in "e"],
            *[f"u_{u}_{i}_m_e" for u in ("1", "3", "4") for i in range(5)],
            *[f"u_{u}_{i}_m_o" for u in ("1", "3", "4") for i in range(3)],
            *[f"u_{u}_{i}_{m}_e" for u in ("2", "6") for i in range(6) for m in ("m", "x")],
            *[f"u_{u}_{i}_x_o" for u in ("2", "6") for i in range(6)],
            *[f"u_{u}_{i}_m_o" for u in ("2", "6") for i in range(3)],
            *[f"u_{u}_{i}_m_e" for u in "5" for i in ("0", "1", "3", "5", "6")],
            *[f"u_{u}_{i}_x_{r}" for u in "5" for i in ("0", "2", "3", "5", "6") for r in ("e", "o")],
            *[f"u_{u}_{i}_m_o" for u in "5" for i in ("0", "2")],
        ]

        energy_params_names = [
            *[rf"$\epsilon_{i}^{{X,{r}}}$" for i in range(2) for r in ("e", "o")],
            *[rf"$\epsilon_{i}^{{M,{r}}}$" for i in "0" for r in ("e", "o")],
            *[rf"$\epsilon_{i}^{{M,{r}}}$" for i in "1" for r in "e"],
            *[rf"$u_{u}^{{{i},e}}$" for u in ("1", "3", "4") for i in range(5)],
            *[rf"$u_{u}^{{{i},o}}$" for u in ("1", "3", "4") for i in range(3)],
            *[rf"$u_{u}^{{{i},{m}e}}$" for u in ("2", "6") for i in range(6) for m in ("M", "X")],
            *[rf"$u_{u}^{{{i},Xo}}$" for u in ("2", "6") for i in range(6)],
            *[rf"$u_{u}^{{{i},Mo}}$" for u in ("2", "6") for i in range(3)],
            *[rf"$u_{u}^{{{i},Me}}$" for u in "5" for i in ("0", "1", "3", "5", "6")],
            *[rf"$u_{u}^{{{i},X{r}}}$" for u in "5" for i in ("0", "2", "3", "5", "6") for r in ("e", "o")],
            *[rf"$u_{u}^{{{i},Mo}}$" for u in "5" for i in ("0", "2")],
        ]

        self._energy_params_dict: Dict[str, FloatParameter] = dict(zip(
            energy_params,
            [FloatParameter(name=p_name) for p_name in energy_params_names]
        ))

        self._general_params_dict: Dict[str, Union[FloatParameter, StringParameter]] = dict(zip(
            ["a", "lamb_m", "lamb_x", "material"],
            [
                FloatParameter(param=1., name=r"$a$"),
                FloatParameter(name=r"$\lambda_M$"),
                FloatParameter(name=r"$\lambda_X$"),
                StringParameter(name="material")
            ]
        ))
        if input_dict is not None:
            self.from_dict(input_dict)

    @property
    def _unique_params_dict(self) -> List[dict]:
        return [self._general_params_dict, self._energy_params_dict]

    @property
    def _protected_params_dict(self) -> Optional[List[dict]]:
        return None

    @property
    def _all_params_dict(self) -> List[dict]:
        return self._unique_params_dict + (self._protected_params_dict or [])

    @property
    def _allowed_params(self):
        return [param_key for param_dict in self._all_params_dict for param_key in param_dict.keys()]

    def _check_key(self, key) -> bool:
        if key not in self._allowed_params:
            warnings.warn(f"Variable {key} is not an expected variable, it is ignored", UserWarning, stacklevel=2)
            return False
        return True

    def __setitem__(self, key, value):
        if self._check_key(key):
            for param_dict in self._unique_params_dict:
                if key in param_dict.keys():
                    param_dict[key].param = value
                    return
            if self._protected_params_dict is not None:
                for param_dict in self._protected_params_dict:
                    if key in param_dict.keys():
                        param_dict[key].param = value
                        warnings.warn(f"The variable {key} is read-only, you should not change it.", UserWarning, stacklevel=2)
                        return
            warnings.warn(f"This should not happen, {key} should be a valid key.", UserWarning, stacklevel=2)

    def __getitem__(self, item) -> Union[float, str]:
        return self.get_param(item) or 0.0

    def get_dict(self) -> dict:
        """Function to get the variables as a dict."""
        out_dict = {}
        for param_dict in self._unique_params_dict:
            for key, item in param_dict.items():
                if item.param is not None:
                    out_dict[key] = item.param
        return out_dict

    def get_param(self, key):
        """Function to get the specific variable"""
        if self._check_key(key):
            for param_dict in self._all_params_dict:
                if key in param_dict.keys():
                    return param_dict[key].param
            warnings.warn(f"This should not happen, {key} should be a valid key.", UserWarning, stacklevel=2)

    def get_name(self, key) -> str:
        """Function to get the name (in LaTeX) of the specific variable"""
        if self._check_key(key):
            for param_dict in self._all_params_dict:
                if key in param_dict.keys():
                    return param_dict[key].name
            warnings.warn(f"This should not happen, {key} should be a valid key.", UserWarning, stacklevel=2)

    def from_dict(self, input_dict: Dict[str, Union[float, str]]):
        """Function to set the variables with a dict.

        Parameters:
            input_dict (dict): The dictionary containing the parameters for the TMD lattice.
                The allowed keys are:
                `a`,            `lamb_m`,       `lamb_x`,       `material`,
                `eps_0_x_e`, `eps_1_x_e`, `eps_0_m_e`, `eps_1_m_e`,
                `u_1_0_m_e`, `u_1_1_m_e`, `u_1_2_m_e`, `u_1_3_m_e`, `u_1_4_m_e`,
                `u_2_0_m_e`, `u_2_1_m_e`, `u_2_2_m_e`, `u_2_3_m_e`, `u_2_4_m_e`, `u_2_5_m_e`,
                `u_2_0_x_e`, `u_2_1_x_e`, `u_2_2_x_e`, `u_2_3_x_e`, `u_2_4_x_e`, `u_2_5_x_e`,
                `u_3_0_m_e`, `u_3_1_m_e`, `u_3_2_m_e`, `u_3_3_m_e`, `u_3_4_m_e`,
                `u_4_0_m_e`, `u_4_1_m_e`, `u_4_2_m_e`, `u_4_3_m_e`, `u_4_4_m_e`,
                `u_5_0_m_e`, `u_5_1_m_e`, `u_5_3_m_e`, `u_5_5_m_e`, `u_5_6_m_e`,
                `u_5_0_x_e`, `u_5_2_x_e`, `u_5_3_x_e`, `u_5_5_x_e`, `u_5_6_x_e`,
                `u_6_0_m_e`, `u_6_1_m_e`, `u_6_2_m_e`, `u_6_3_m_e`, `u_6_4_m_e`, `u_6_5_m_e`,
                `u_6_0_x_e`, `u_6_1_x_e`, `u_6_2_x_e`, `u_6_3_x_e`, `u_6_4_x_e`, `u_6_5_x_e`,
                `eps_0_x_o`, `eps_1_x_o`, `eps_0_m_o`,
                `u_1_0_m_o`, `u_1_1_m_o`, `u_1_2_m_o`,
                `u_2_0_m_o`, `u_2_1_m_o`, `u_2_2_m_o`,
                `u_2_0_x_o`, `u_2_1_x_o`, `u_2_2_x_o`, `u_2_3_x_o`, `u_2_4_x_o`, `u_2_5_x_o`,
                `u_3_0_m_o`, `u_3_1_m_o`, `u_3_2_m_o`,
                `u_4_0_m_o`, `u_4_1_m_o`, `u_4_2_m_o`,
                `u_5_0_m_o`, `u_5_2_m_o`,
                `u_5_0_x_o`, `u_5_2_x_o`, `u_5_3_x_o`, `u_5_5_x_o`, `u_5_6_x_o`,
                `u_6_0_m_o`, `u_6_1_m_o`, `u_6_2_m_o`,
                `u_6_0_x_o`, `u_6_1_x_o`, `u_6_2_x_o`, `u_6_3_x_o`, `u_6_4_x_o` and `u_6_5_x_o`
        """
        for param_name, value in input_dict.items():
            self[param_name] = value


_liu_2nn_mos2_fitted = ParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.318955279,    0.073,          None,           "MoS2",
        -4.75202031,    -3.81191516,
        -0.18254965,    0.56039034,     -0.35035407,    0.02587717,     0.32520158,     0.22190189
    ]
)))

_liu_6nn_mos2_fitted = ParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.318955279,    0.073,          None,           "MoS2",
        -5.09750973,    -4.10128918,
        -0.14342346,    0.50944507,     0.11371639,     0.07975301,     0.16337675,     0.08534869,
        0.05829675,     -0.07413799,    -0.04025448,    0.17963519,     0.26548849,
        -0.03801049,    0.00373271,     -0.04490756,    -0.15506986,    -0.17736107,    0.26986353
    ]
)))

_liu_2nn_mos2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.319,          0.073,          None,           "MoS2",
        1.046,          2.104,
        -0.184,         0.507,          -0.401,         0.057,          0.338,          0.218
    ]
)))

_liu_2nn_mose2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3326,         0.091,          None,           "MoSe2",
        0.919,          2.065,
        -0.188,         0.456,          -0.317,         0.13,           0.29,           0.211
    ]
)))

_liu_2nn_mote2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3557,         0.107,          None,           "MoTe2",
        0.605,          1.972,
        -0.169,         0.39,           -0.228,         0.252,          0.239,          0.207,
    ]
)))

_liu_2nn_ws2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3191,         0.211,          None,           "WS2",
        1.13,           2.275,
        -0.206,         0.536,          -0.567,        -0.061,         0.384,          0.286
    ]
)))

_liu_2nn_wse2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3325,         0.228,          None,           "WSe2",
        0.943,          2.179,
        -0.207,         0.486,          -0.457,         0.034,         0.329,            0.263
    ]
)))

_liu_2nn_wte2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.356,          0.237,          None,           "WTe2",
        0.606,          2.102,
        -0.175,         0.41,           -0.342,         0.19,           0.27,           0.233
    ]
)))

_liu_6nn_mos2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3190,         0.073,          None,           "MoS2",
        0.683,          1.707,
        -0.146,         0.506,          0.114,          0.073,          0.162,          0.085,
        0.060,          -0.077,         -0.034,         0.167,          0.273,
        -0.038,         0.001,          -0.046,         -0.150,         -0.176,         0.266
    ]
)))

_liu_6nn_mose2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3326,         0.091,          None,           "MoSe2",
        0.684,          1.546,
        -0.146,         0.432,          0.13,           0.075,          0.117,          0.144,
        0.039,          -0.0797,        0.0174,         0.1559,         0.2413,
        -0.042,         0.008,          -0.036,         -0.15,          -0.172,         0.272,
    ]
)))

_liu_6nn_mote2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3557,         0.107,          None,           "MoTe2",
        0.588,          1.303,
        -0.226,         0.036,          0.234,          0.017,          0.098,          0.4,
        0.003,          0.1951,         0.0526,         0.1703,         0.0289,
        0.057,          0.187,          -0.103,         0.087,          -0.141,         -0.045,
    ]
)))

_liu_6nn_ws2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3191,         0.211,          None,           "WS2",
        0.717,          1.916,
        -0.152,         0.59,           0.097,          0.016,           0.178,          0.047,
        0.069,          -0.1236,        -0.0659,        0.1858,          0.3014,
        -0.054,         0.002,          -0.045,         -0.163,          -0.206,         0.325,
    ]
)))

_liu_6nn_wse2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3325,         0.228,          None,           "WSe2",
        0.728,          1.655,
        -0.146,         0.507,          0.124,          0.015,          0.127,          0.117,
        0.036,          -0.1236,        0.0007,         0.1739,         0.2702,
        -0.061,         0.007,          -0.032,         -0.164,         -0.202,         0.329,
    ]
)))

_liu_6nn_wte2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.356,          0.237,          None,           "WTe2",
        0.697,          1.38,
        -0.109,         0.368,          0.164,          0.038,          0.093,           0.204,
        -0.015,         -0.1236,        0.110,          0.1306,         0.241,
        -0.066,         -0.013,         -0.011,        -0.132,          -0.177,         0.312
    ]
)))

_wu = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.085433 and https://doi.org/10.1103/PhysRevB.91.075310
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
        -0.146,         0.506,          0.114,         0.073,          0.162,         0.085,
        -0.189,         -0.024,         -0.117,
        0.060,          -0.077,         -0.034,         0.167,          0.273,
        -0.063,         0.025,
        -0.038,         0.001,          -0.046,          -0.150,        -0.176,          0.266,
        0.165,          0.140,          -0.122
    ]
)))

_fang_mos2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.318,          0.0836,         0.0556,         "MoS2",
        -0.138,         0.0874,         1.0688,
        -1.9065,        -2.8949,        -1.2902,        -0.7755,
        1.4114,         -0.9402,        0.6517,         -0.8836,        -0.9535,
        -0.7883,        2.1584,         -1.379,
        -0.2979,        0.4096,         -0.1145,        -0.5581,        0.2487,         0.2747,
        -0.2069,        -0.2562,        0.0323,
        0.9122,         -0.0385,        -0.1063,        0.0059,         0.0075,         -0.1916,
        0.8651,         -0.0705,        0.0995,         -0.1872,        -0.0679,         -0.1739,
        None,           -0.1498,        -0.2451,        -0.0686,        -0.2205
    ]
)))

_fang_ws2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.318,          0.2874,         0.0556,         "WS2",
        -0.0393,        0.1984,         1.3754,
        -2.3461,        -3.3706,        -1.5534,        -1.1278,
        1.5629,         -0.9878,        0.6718,         -1.013,         -0.9491,
        -0.8855,        2.3121,         -1.4376,
        -0.3716,        0.4896,         -0.1467,        -0.6892,         0.303,         0.3537,
        -0.2011,        -0.3106,        0.0263,
        0.9673,         -0.1018,        -0.1645,         0.0143,         -0.0315,         -0.2112,
        0.8726,         -0.0989,        0.1105,         -0.2187,        -0.0818,         -0.1749,
        None,           -0.1533,        -0.2736,        -0.0659,        -0.2618
    ]
)))

_fang_mose2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.332,          0.0836,         0.247,          "MoSe2",
        -0.2297,        0.0149,         0.7819,
        -1.7806,        -2.9015,        -1.1726,        -0.6567,
        1.2677,         -0.8738,        0.5545,         -0.772,         -0.8578,
        -0.6946,        1.9415,         -1.3258,
        -0.2636,        0.352,          -0.096,         -0.4734,        0.2012,         0.2505,
        -0.146,         -0.1912,        0.0177,
        0.9911,         -0.0394,        -0.1216,        -0.0036,        0.0047,         -0.2166,
        0.9638,         -0.068,         0.0755,         -0.1724,        -0.0735,        -0.2112,
        None,           -0.1553,        -0.2154,        -0.0691,        -0.2227
    ]
)))

_fang_wse2 = ParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.332,          0.2874,         0.247,          "WSe2",
        -0.1667,        0.0984,         1.0349,
        -2.182,         -3.3642,        -1.3937,        -0.9573,
        1.403,          -0.9044,        0.5711,         -0.8998,        -0.8548,
        -0.7744,        2.0858,         -1.4014,
        -0.333,         0.4233,         -0.125,          -0.5837,        0.2456,        0.319,
        -0.1395,        -0.2321,         0.0129,
        1.047,          -0.1027,         -0.1857,         0.0029,         -0.0377,        -0.2399,
        0.9763,         -0.092,          0.0797,        -0.1985,        -0.0912,        -0.2171,
        None,           -0.1608,        -0.2424,        -0.0676,        -0.2618
    ]
)))


_wu_fitted = ParametersList(dict(zip(
    #  fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
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
        -5.098,         -4.101,         -2.246,
        -0.143,         0.509,           0.114,         0.079,           0.164,         0.085,
        -0.167,         -0.022,         -0.148,
        0.058,          -0.074,         -0.040,         0.179,          0.265,
        -0.061,         0.037,
        -0.038,         0.004,          -0.045,         -0.155,         -0.177,          0.270,
        0.169,          0.149,          -0.123
    ]
)))

_fang_mos2_fitted = ParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.3189552789,    0.0836,         0.0556,         "MoS2",
        -6.48602262,    -5.18473506,    -4.59364376,
        -7.7577787,     -9.05082632,    -7.23255442,    -7.02669394,
        1.26825943,     -1.08658995,    0.73122123,     -0.69014034,    -0.84531854,
        -0.78156314,    2.19543484,     -1.3171294,
        -0.03516888,    0.47269242,     -0.12368948,     -0.39732265,    0.25919776,    0.14680296,
        -0.16375119,    0.10644792,    -0.08866745,
        0.87185047,     0.13267939,     -0.06673336,     -0.08669893,    -0.1665462,     -0.25628203,
        0.77863574,     0.05642279,     0.055016,      -0.07622647,    0.01287608,     -0.09818404,
        None,           -0.30871332,    -0.14419783,    0.01966645,     -0.37086515
    ]
)))

_jorissen_mos2 = ParametersList(dict(zip(
    # https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "eps_0_x_e",    "eps_1_x_e",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",

    ],
    [
        0.3189552789,    0.0836,         0.0556,         "MoS2",
        -6.47506922,    -4.89138157,
        -7.90697285,    -9.47021899,
        0.99874101,     -1.28945959,    0.79526257,     -0.68815634,    -0.79456582,
        -0.04758008,    0.58038247,     0.07374717,    -0.41375629,      0.29861063,    0.04452479,
        0.79492166,     0.24849791,    -0.16415515,     -0.0019305,     -0.2934992,     -0.17447978,
    ]
)))

_all_mos2 = ParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_x_e", "eps_1_x_e", "eps_0_m_e", "eps_1_m_e",
        "u_1_0_m_e", "u_1_1_m_e", "u_1_2_m_e", "u_1_3_m_e", "u_1_4_m_e",
        "u_2_0_m_e", "u_2_1_m_e", "u_2_2_m_e", "u_2_3_m_e", "u_2_4_m_e", "u_2_5_m_e",
        "u_2_0_x_e", "u_2_1_x_e", "u_2_2_x_e", "u_2_3_x_e", "u_2_4_x_e", "u_2_5_x_e",
        "u_3_0_m_e", "u_3_1_m_e", "u_3_2_m_e", "u_3_3_m_e", "u_3_4_m_e",
        "u_4_0_m_e", "u_4_1_m_e", "u_4_2_m_e", "u_4_3_m_e", "u_4_4_m_e",
        "u_5_0_m_e", "u_5_1_m_e", "u_5_3_m_e", "u_5_5_m_e", "u_5_6_m_e",
        "u_5_0_x_e", "u_5_2_x_e", "u_5_3_x_e", "u_5_5_x_e", "u_5_6_x_e",
        "u_6_0_m_e", "u_6_1_m_e", "u_6_2_m_e", "u_6_3_m_e", "u_6_4_m_e", "u_6_5_m_e",
        "u_6_0_x_e", "u_6_1_x_e", "u_6_2_x_e", "u_6_3_x_e", "u_6_4_x_e", "u_6_5_x_e",
        "eps_0_x_o", "eps_1_x_o", "eps_0_m_o",
        "u_1_0_m_o", "u_1_1_m_o", "u_1_2_m_o",
        "u_2_0_m_o", "u_2_1_m_o", "u_2_2_m_o",
        "u_2_0_x_o", "u_2_1_x_o", "u_2_2_x_o", "u_2_3_x_o", "u_2_4_x_o", "u_2_5_x_o",
        "u_3_0_m_o", "u_3_1_m_o", "u_3_2_m_o",
        "u_4_0_m_o", "u_4_1_m_o", "u_4_2_m_o",
        "u_5_0_m_o", "u_5_2_m_o",
        "u_5_0_x_o", "u_5_2_x_o", "u_5_3_x_o", "u_5_5_x_o", "u_5_6_x_o",
        "u_6_0_m_o", "u_6_1_m_o", "u_6_2_m_o",
        "u_6_0_x_o", "u_6_1_x_o", "u_6_2_x_o", "u_6_3_x_o", "u_6_4_x_o", "u_6_5_x_o",

    ],
    [
        0.3189552789,    0.0836,         0.0556,         "MoS2",
        -5.85535089e+00, -6.55002800e+00, -4.11031485e+00, -3.50757649e+00,
        1.35813632e+00, -9.39814832e-01, 6.36315099e-01, -8.65807106e-01, -9.40867659e-01,
        -2.88698499e-01, 5.08367498e-01, -1.17988537e-01, -5.53682676e-01, 2.28504870e-01, 2.29751185e-01,
        9.10803109e-01, 2.24402452e-02, -9.61186364e-02, 3.12400382e-03, -2.79471707e-02, -1.84215453e-01,
        8.19856098e-02, -1.36942029e-01, -2.32226175e-01, -7.56614354e-02, -2.18453245e-01,
        -1.05966154e-02, 8.18751691e-04, -3.19735380e-02, 1.26355104e-03, 2.88724090e-02,
        5.40497529e-04, -7.74247515e-03, 3.39424467e-02, 9.11337396e-03, -5.23766300e-03,
        -1.31303857e-02, -2.16101527e-05, 2.18892380e-03, 2.85609057e-03, -5.35732464e-03,
        2.67882896e-04, -2.91376849e-02, -1.15194554e-04, -5.58489128e-02, -6.08538922e-02, 7.25848009e-02,
        4.85565562e-02, -6.95112283e-02, -6.91879264e-03, 3.04557146e-02, 6.86382071e-02, -1.97025368e-03,
        -4.96490884e+00, -4.73278639e+00, -2.81154976e+00,
        -7.82104656e-01, 2.13753469e+00, -1.45504793e+00,
        -2.62365349e-01, -2.15615929e-01, 1.02434139e-02,
        8.49785326e-01, -7.06226486e-02, 4.13053364e-02, -1.12769330e-01, 4.11554466e-02, -1.73215911e-01,
        3.14319913e-02, -6.35638918e-02, 5.80612006e-02,
        -3.35710022e-02, 6.33229092e-02, 5.05315655e-02,
        -6.50386867e-02, 1.28907505e-01,
        -2.72518553e-02, -8.45574502e-03, -7.19418647e-02, 9.94491714e-03, -2.39873085e-02,
        4.30630790e-02, 5.98009509e-02, -1.26888729e-01,
        7.30968724e-02, -4.37471809e-02, -9.40970153e-03, 2.87364009e-02, -2.88204369e-02, 2.28460254e-02
    ]
)))

liu2 = {
    "MoS2": _liu_2nn_mos2,
    "MoSe2": _liu_2nn_mose2,
    "MoTe2": _liu_2nn_mote2,
    "WS2": _liu_2nn_ws2,
    "WSe2": _liu_2nn_wse2,
    "WTe2": _liu_2nn_wte2,
    "fitted": {"MoS2": _liu_2nn_mos2_fitted}
}

liu6 = {
    "MoS2": _liu_6nn_mos2,
    "MoSe2": _liu_6nn_mose2,
    "MoTe2": _liu_6nn_mote2,
    "WS2": _liu_6nn_ws2,
    "WSe2": _liu_6nn_wse2,
    "WTe2": _liu_6nn_wte2,
    "fitted": {"MoS2": _liu_6nn_mos2_fitted}
}

wu = {
    "MoS2": _wu,
    "fitted": {"MoS2": _wu_fitted}
}

fang = {
    "MoS2": _fang_mos2,
    "MoSe2": _fang_mose2,
    "WS2": _fang_ws2,
    "WSe2": _fang_wse2,
    "fitted": {"MoS2": _fang_mos2_fitted}
}

jorissen = {
    "MoS2": _jorissen_mos2
}

all = {
    "MoS2": _all_mos2
}
