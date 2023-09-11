"""Class to save the parameters for TMD models"""
from dataclasses import dataclass, fields
from typing import Optional, Dict, Union, List
import warnings
import numpy as np

@dataclass
class Parameter:
    """Class to store one seperate parameter"""
    name: str = ""

@dataclass
class FloatParameter(Parameter):
    param: Optional[float] = None

@dataclass
class StringParameter(Parameter):
    param: Optional[str] = None


class ParametersList:
    """Class to save the parameters"""

    def __init__(self, input_dict: Optional[Dict[str, Union[float, str]]] = None):

        self._general_params = ["a", "lamb_m", "lamb_x", "material"]
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
            *[rf"$u_{u}_{{{i},Mo}}$" for u in "5" for i in ("0", "2")],
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
    def _allowed_params(self):
        return [*self._general_params_dict, *self._energy_params_dict]

    def _check_key(self, key) -> bool:
        if not key in self._allowed_params:
            warnings.warn("Variable {0} is not an expected variable, it is ignored".format(
                key
            ), UserWarning, stacklevel=2)
            return False
        return True
    def __setitem__(self, key, value):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                self._general_params_dict[key].param = value
            elif key in [*self._energy_params_dict]:
                self._energy_params_dict[key].param = value
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def __getitem__(self, item) -> Union[float, str]:
        return self.get_param(item) or 0.0

    def get_param(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].param
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].param
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def get_name(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].name
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].name
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def from_dict(self, input_dict: Dict[str, Union[float, str]]):
        """Function to set the variables with a dict"""
        for param_name, value in input_dict.items():
            self[param_name] = value


class SKParametersList(ParametersList):
    """Class to store the parameters for SK e/o model"""

    def __init__(self, input_dict: Optional[Dict[str, Union[str, Optional[float]]]] = None):
        super().__init__(None)
        self._recalculate_params_bool: bool = True
        self._tan_theta: Optional[float] = None
        sk_params = [
            "theta",
            *[f"delta_{r}" for r in ("p", "z", "0", "1", "2")],
            *[f"v_0_pp{r}" for r in ("s", "p")],
            *[f"v_{n}_{r}_pd{i}" for n in ("1", "3", "4") for r in ("e", "o") for i in ("s", "p")],
            *[f"v_{n}_{r}_pp{i}{t}" for n in ("2", "5", "6") for r in ("e", "o") for i in ("s", "p") for t in ("", "_tb")],
            *[f"v_{n}_e_dds" for n in ("2", "5", "6")],
            *[f"v_{n}_{r}_dd{i}" for n in ("2", "5", "6") for r in ("e", "o") for i in ("p", "d")]
        ]

        sk_params_names = [
            r"$\theta$",
            *[rf"$\delta_{r}$" for r in ("p", "z", "0", "1", "2")],
            *[rf"$V^0_{{pp{r}}}$" for r in (r"\sigma", r"\pi")],
            *[rf"$V^{{{n}_{r}}}_{{pd{i}}}$" for n in ("1", "3", "4") for r in ("e", "o") for i in (r"\sigma", r"\pi")],
            *[rf"$V^{{{n}_{r}}}_{{pp{i}{t}}}$" for n in ("2", "5", "6") for r in ("e", "o") for i in (r"\sigma", r"\pi") for t in ("", ",tb")],
            *[rf"$V^{{{n}_e}}_{{dd\sigmas}}$" for n in ("2", "5", "6")],
            *[rf"$V^{{{n}_{r}}}_{{dd{i}}}$" for n in ("2", "5", "6") for r in ("e", "o") for i in (r"\pi", r"\delta")]
        ]
        self._sk_params_dict: dict = dict(zip(
            sk_params,
            [FloatParameter(name=p_name) for p_name in sk_params_names]
        ))
        if input_dict is not None:
            self.from_dict(input_dict)

    def from_dict(self, input_dict: Dict[str, Union[float, str]]):
        self._recalculate_params_bool = False
        super().from_dict(input_dict)
        self._recalculate_params_bool = True
        self._recalculate_params()

    def _recalculate_params(self):
        if self._recalculate_params_bool:
            self._tan_theta = np.tan(self["theta"]) if self["theta"] is not None else np.sqrt(3 / 4)
            # Onsite elements
            # X
            self._energy_params_dict["eps_0_x_e"].param = self._comb_nonefloat([self["delta_p"], self["v_0_ppp"]])
            self._energy_params_dict["eps_1_x_e"].param = self._subtract_param(self["delta_z"], self["v_0_pps"])
            self._energy_params_dict["eps_0_x_o"].param = self._subtract_param(self["delta_p"], self["v_0_ppp"])
            self._energy_params_dict["eps_1_x_o"].param = self._comb_nonefloat([self["delta_z"], self["v_0_pps"]])
            # M
            self._energy_params_dict["eps_0_m_e"].param = self["delta_0"]
            self._energy_params_dict["eps_1_m_e"].param = self["delta_2"]
            self._energy_params_dict["eps_0_m_o"].param = self["delta_1"]
            # h1
            [
                self._energy_params_dict["u_1_0_m_e"].param,
                self._energy_params_dict["u_1_1_m_e"].param,
                self._energy_params_dict["u_1_2_m_e"].param,
                self._energy_params_dict["u_1_3_m_e"].param,
                self._energy_params_dict["u_1_4_m_e"].param
            ] = self._h_mx_e(r=-1, vpds=self["v_1_e_pds"], vpdp=self["v_1_e_pdp"])
            [
                self._energy_params_dict["u_1_0_m_o"].param,
                self._energy_params_dict["u_1_1_m_o"].param,
                self._energy_params_dict["u_1_2_m_o"].param
            ] = self._h_mx_o(r=-1, vpds=self["v_1_o_pds"], vpdp=self["v_1_o_pdp"])
            # h2
            [
                self._energy_params_dict["u_2_0_m_e"].param,
                self._energy_params_dict["u_2_1_m_e"].param,
                self._energy_params_dict["u_2_2_m_e"].param,
                self._energy_params_dict["u_2_3_m_e"].param,
                self._energy_params_dict["u_2_4_m_e"].param,
                self._energy_params_dict["u_2_5_m_e"].param
            ] = self._h_mm_x_e(vdds=self["v_2_e_dds"], vddp=self["v_2_e_ddp"], vddd=self["v_2_e_ddd"])
            [
                self._energy_params_dict["u_2_0_m_o"].param,
                self._energy_params_dict["u_2_1_m_o"].param,
                self._energy_params_dict["u_2_2_m_o"].param
            ] = self._h_mm_x_o(vddp=self["v_2_o_ddp"], vddd=self["v_2_o_ddd"])
            [
                self._energy_params_dict["u_2_0_x_e"].param,
                self._energy_params_dict["u_2_1_x_e"].param,
                self._energy_params_dict["u_2_2_x_e"].param,
                self._energy_params_dict["u_2_3_x_e"].param,
                self._energy_params_dict["u_2_4_x_e"].param,
                self._energy_params_dict["u_2_5_x_e"].param
            ] = self._h_xx_x_e(r=np.sqrt(3), vpps=self["v_2_e_pps"], vppp=self["v_2_e_ppp"],
                               vppstb=self["v_2_e_pps_tb"], vppptb=self["v_2_e_ppp_tb"])
            [
                self._energy_params_dict["u_2_0_x_o"].param,
                self._energy_params_dict["u_2_1_x_o"].param,
                self._energy_params_dict["u_2_2_x_o"].param,
                self._energy_params_dict["u_2_3_x_o"].param,
                self._energy_params_dict["u_2_4_x_o"].param,
                self._energy_params_dict["u_2_5_x_o"].param
            ] = self._h_xx_x_o(r=np.sqrt(3), vpps=self["v_2_o_pps"], vppp=self["v_2_o_ppp"],
                               vppstb=self["v_2_o_pps_tb"], vppptb=self["v_2_o_ppp_tb"])
            # h3
            [
                self._energy_params_dict["u_3_0_m_e"].param,
                self._energy_params_dict["u_3_1_m_e"].param,
                self._energy_params_dict["u_3_2_m_e"].param,
                self._energy_params_dict["u_3_3_m_e"].param,
                self._energy_params_dict["u_3_4_m_e"].param
            ] = self._h_mx_e(r=2, vpds=self["v_3_e_pds"], vpdp=self["v_3_e_pdp"])
            [
                self._energy_params_dict["u_3_0_m_o"].param,
                self._energy_params_dict["u_3_1_m_o"].param,
                self._energy_params_dict["u_3_2_m_o"].param
            ] = self._h_mx_o(r=2, vpds=self["v_3_o_pds"], vpdp=self["v_3_o_pdp"])
            # h4
            [
                self._energy_params_dict["u_4_0_m_e"].param,
                self._energy_params_dict["u_4_1_m_e"].param,
                self._energy_params_dict["u_4_2_m_e"].param,
                self._energy_params_dict["u_4_3_m_e"].param,
                self._energy_params_dict["u_4_4_m_e"].param
            ] = self._h_mx_e(r=-np.sqrt(7), vpds=self["v_4_e_pds"], vpdp=self["v_4_e_pdp"])
            [
                self._energy_params_dict["u_4_0_m_o"].param,
                self._energy_params_dict["u_4_1_m_o"].param,
                self._energy_params_dict["u_4_2_m_o"].param
            ] = self._h_mx_o(r=-np.sqrt(7), vpds=self["v_4_o_pds"], vpdp=self["v_4_o_pdp"])
            # h5
            [
                self._energy_params_dict["u_5_0_m_e"].param,
                self._energy_params_dict["u_5_1_m_e"].param,
                self._energy_params_dict["u_5_3_m_e"].param,
                self._energy_params_dict["u_5_5_m_e"].param,
                self._energy_params_dict["u_5_6_m_e"].param
            ] = self._h_mm_y_e(vdds=self["v_5_e_dds"], vddp=self["v_5_e_ddp"], vddd=self["v_5_e_ddd"])
            [
                self._energy_params_dict["u_5_0_m_o"].param,
                self._energy_params_dict["u_5_2_m_o"].param
            ] = self._h_mm_y_o(vddp=self["v_5_o_ddp"], vddd=self["v_5_o_ddd"])
            [
                self._energy_params_dict["u_5_0_x_e"].param,
                self._energy_params_dict["u_5_2_x_e"].param,
                self._energy_params_dict["u_5_3_x_e"].param,
                self._energy_params_dict["u_5_5_x_e"].param,
                self._energy_params_dict["u_5_6_x_e"].param
            ] = self._h_xx_y_e(r=3, vpps=self["v_5_e_pps"], vppp=self["v_5_e_ppp"],
                               vppstb=self["v_5_e_pps_tb"], vppptb=self["v_5_e_ppp_tb"])
            [
                self._energy_params_dict["u_5_0_x_o"].param,
                self._energy_params_dict["u_5_2_x_o"].param,
                self._energy_params_dict["u_5_3_x_o"].param,
                self._energy_params_dict["u_5_5_x_o"].param,
                self._energy_params_dict["u_5_6_x_o"].param
            ] = self._h_xx_y_o(r=3, vpps=self["v_5_o_pps"], vppp=self["v_5_o_ppp"],
                               vppstb=self["v_5_e_pps_tb"], vppptb=self["v_5_o_ppp_tb"])
            # h6
            [
                self._energy_params_dict["u_6_0_m_e"].param,
                self._energy_params_dict["u_6_1_m_e"].param,
                self._energy_params_dict["u_6_2_m_e"].param,
                self._energy_params_dict["u_6_3_m_e"].param,
                self._energy_params_dict["u_6_4_m_e"].param,
                self._energy_params_dict["u_6_5_m_e"].param
            ] = self._h_mm_x_e(vdds=self["v_6_e_dds"], vddp=self["v_6_e_ddp"], vddd=self["v_6_e_ddd"])
            [
                self._energy_params_dict["u_6_0_m_o"].param,
                self._energy_params_dict["u_6_1_m_o"].param,
                self._energy_params_dict["u_6_2_m_o"].param
            ] = self._h_mm_x_o(vddp=self["v_6_o_ddp"], vddd=self["v_6_o_ddd"])
            [
                self._energy_params_dict["u_6_0_x_e"].param,
                self._energy_params_dict["u_6_1_x_e"].param,
                self._energy_params_dict["u_6_2_x_e"].param,
                self._energy_params_dict["u_6_3_x_e"].param,
                self._energy_params_dict["u_6_4_x_e"].param,
                self._energy_params_dict["u_6_5_x_e"].param
            ] = self._h_xx_x_e(r=np.sqrt(3), vpps=self["v_6_e_pps"], vppp=self["v_6_e_ppp"],
                               vppstb=self["v_6_e_pps_tb"], vppptb=self["v_6_e_ppp_tb"])
            [
                self._energy_params_dict["u_6_0_x_o"].param,
                self._energy_params_dict["u_6_1_x_o"].param,
                self._energy_params_dict["u_6_2_x_o"].param,
                self._energy_params_dict["u_6_3_x_o"].param,
                self._energy_params_dict["u_6_4_x_o"].param,
                self._energy_params_dict["u_6_5_x_o"].param
            ] = self._h_xx_x_o(r=np.sqrt(3), vpps=self["v_6_o_pps"], vppp=self["v_6_o_ppp"],
                               vppstb=self["v_6_o_pps_tb"], vppptb=self["v_6_o_ppp_tb"])


    def _subtract_param(self, p1: Optional[float], p2: Optional[float]):
        return self._comb_nonefloat([p1, -p2 if p2 is not None else None])

    @staticmethod
    def _h_mm_x_e(vdds: Optional[float], vddp: Optional[float], vddd: Optional[float]) -> List[Optional[float]]:
        vddsd_bool = (vdds is not None) and (vddd is not None)
        vddp_bool = vddp is not None
        u_0: Optional[float] = 1 / 4 * (vdds + 3 * vddd) if vddsd_bool else None
        u_1: Optional[float] = np.sqrt(3) / 4 * (-vdds + vddd) if vddsd_bool else None
        u_2: Optional[float] = 0 if (vddsd_bool or vddp_bool) else None
        u_3: Optional[float] = 1 / 4 * (3 * vdds + vddd) if vddsd_bool else None
        u_4: Optional[float] = 0 if (vddsd_bool or vddp_bool) else None
        u_5: Optional[float] = vddp if vddp_bool else None
        return [u_0, u_1, u_2, u_3, u_4, u_5]

    @staticmethod
    def _h_mm_y_e(vdds: Optional[float], vddp: Optional[float], vddd: Optional[float]) -> List[Optional[float]]:
        vddsd_bool = (vdds is not None) and (vddd is not None)
        vddp_bool = vddp is not None
        u_0: Optional[float] = 1 / 4 * (vdds + 3 * vddd) if vddsd_bool else None
        u_1: Optional[float] = np.sqrt(3) / 4 * (-vdds + vddd) if vddsd_bool else None
        u_3: Optional[float] = 1 / 4 * (3 * vdds + vddd) if vddsd_bool else None
        u_5: Optional[float] = vddp if vddp_bool else None
        u_6: Optional[float] = np.sqrt(3) / 4 * (-vdds + vddd) if vddsd_bool else None
        return [u_0, u_1, u_3, u_5, u_6]

    @staticmethod
    def _h_mm_x_o(vddp: Optional[float], vddd: Optional[float]) -> List[Optional[float]]:
        vddd_bool = vddd is not None
        vddp_bool = vddp is not None
        u_0: Optional[float] = vddp if vddp_bool else None
        u_1: Optional[float] = 0 if (vddd_bool or vddp_bool) else None
        u_2: Optional[float] = vddd if vddd_bool else None
        return [u_0, u_1, u_2]

    @staticmethod
    def _h_mm_y_o(vddp: Optional[float], vddd: Optional[float]) -> List[Optional[float]]:
        vddd_bool = vddd is not None
        vddp_bool = vddp is not None
        u_0: Optional[float] = vddp if vddp_bool else None
        u_2: Optional[float] = vddd if vddd_bool else None
        return [u_0, u_2]

    @staticmethod
    def _comb_nonefloat(u_list: List[Optional[float]]) -> Optional[float]:
        bl = [ul is not None for ul in u_list]
        return None if np.all(np.logical_not(bl)) else np.sum([u_list[li] for li in np.arange(len(u_list))[bl]])

    def _h_xx_x_e(self, r: float, vpps: Optional[float], vppp: Optional[float],
                  vppstb: Optional[float] = None, vppptb: Optional[float] = None) -> List[Optional[float]]:
        vpps_bool = vpps is not None
        vppp_bool = vppp is not None
        vptb_bool = vppstb is not None
        vstb_bool = vppptb is not None
        u_0_0: Optional[float] = vpps if vpps_bool else None
        u_0_1: Optional[float] = vppptb if vptb_bool else None
        u_0_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_1: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_4: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = -vppstb if vstb_bool else None
        u_5_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        return [u_0, u_1, u_2, u_3, u_4, u_5]

    def _h_xx_y_e(self, r: float, vpps: Optional[float], vppp: Optional[float],
                  vppstb: Optional[float] = None, vppptb: Optional[float] = None) -> List[Optional[float]]:
        vpps_bool = vpps is not None
        vppp_bool = vppp is not None
        vptb_bool = vppstb is not None
        vstb_bool = vppptb is not None
        u_0_0: Optional[float] = vpps if vpps_bool else None
        u_0_1: Optional[float] = vppptb if vptb_bool else None
        u_0_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = -vppstb if vstb_bool else None
        u_5_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        u_6: Optional[float] = 2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        return [u_0, u_2, u_3, u_5, u_6]

    def _h_xx_x_o(self, r: float, vpps: Optional[float], vppp: Optional[float],
                  vppstb: Optional[float] = None, vppptb: Optional[float] = None) -> List[Optional[float]]:
        vpps_bool = vpps is not None
        vppp_bool = vppp is not None
        vptb_bool = vppstb is not None
        vstb_bool = vppptb is not None
        u_0_0: Optional[float] = vpps if vpps_bool else None
        u_0_1: Optional[float] = -vppptb if vptb_bool else None
        u_0_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_1: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = -vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_4: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = vppstb if vstb_bool else None
        u_5_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        return [u_0, u_1, u_2, u_3, u_4, u_5]

    def _h_xx_y_o(self, r: float, vpps: Optional[float], vppp: Optional[float],
                  vppstb: Optional[float] = None, vppptb: Optional[float] = None) -> List[Optional[float]]:
        vpps_bool = vpps is not None
        vppp_bool = vppp is not None
        vptb_bool = vppstb is not None
        vstb_bool = vppptb is not None
        u_0_0: Optional[float] = vpps if vpps_bool else None
        u_0_1: Optional[float] = -vppptb if vptb_bool else None
        u_0_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_2: Optional[float] = 2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = -vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = vppstb if vstb_bool else None
        u_5_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        u_6: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2) if vstb_bool and vptb_bool else None
        return [u_0, u_2, u_3, u_5, u_6]

    def _h_mx_e(self, r: float, vpds: Optional[float], vpdp: Optional[float]) -> List[Optional[float]]:
        if (vpds is not None) and (vpdp is not None):
            u_0 = np.sqrt(2) * r * vpdp / (np.sqrt(r ** 2 + self._tan_theta ** 2))
            u_1 = -r * (2 * np.sqrt(3) * self._tan_theta ** 2 * vpdp + (r ** 2 - 2 * self._tan_theta ** 2) * vpds) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_2 = -r * (2 * self._tan_theta ** 2 * vpdp + np.sqrt(3) * r ** 2 * vpds) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_3 = self._tan_theta * (r ** 2 * 2 * np.sqrt(3) * vpdp + (2 * self._tan_theta ** 2 - r ** 2) * vpds) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_4 = r ** 2 * self._tan_theta * (2 * vpdp - np.sqrt(3) * vpds) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            return [u_0, u_1, u_2, u_3, u_4]
        else:
            return [None, None, None, None, None]

    def _h_mx_o(self, r: float, vpds: Optional[float], vpdp: Optional[float]) -> List[Optional[float]]:
        if (vpds is not None) and (vpdp is not None):
            u_0 = np.sqrt(2) * self._tan_theta * vpdp / np.sqrt(r ** 2 + self._tan_theta ** 2)
            u_1 = np.sqrt(2) * self._tan_theta * ((self._tan_theta ** 2 - r ** 2) * vpdp + r ** 2 * np.sqrt(3) * vpds) / ((r ** 2 + self._tan_theta ** 2) ** (3 / 2))
            u_2 = np.sqrt(2) * r * ((r ** 2 - self._tan_theta ** 2) * vpdp + self._tan_theta ** 2 * np.sqrt(3) * vpds) / ((r ** 2 + self._tan_theta ** 2) ** (3 / 2))
            return [u_0, u_1, u_2]
        else:
            return [None, None, None]


    @property
    def _allowed_params(self):
        return [*self._general_params_dict, *self._energy_params_dict, *self._sk_params_dict]

    def __setitem__(self, key, value):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                self._general_params_dict[key].param = value
            elif key in [*self._sk_params_dict]:
                self._sk_params_dict[key].param = value
                self._recalculate_params()
            elif key in [*self._energy_params_dict]:
                warnings.warn("The variable {0} is set by SK functions and is read-only, don't change it.".format(key))
                self._energy_params_dict[key].param = value
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def get_param(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].param
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].param
            elif key in [*self._sk_params_dict]:
                return self._sk_params_dict[key].param
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def get_name(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].name
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].name
            elif key in [*self._sk_params_dict]:
                return self._sk_params_dict[key].name
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))


class SKSimpleParametersList(SKParametersList):
    def __init__(self, input_dict: Optional[Dict[str, Union[str, Optional[float]]]] = None):
        super().__init__(None)
        sk_simple_params = [
            *[f"v_{n}_pd{i}" for n in ("1", "3", "4") for i in ("s", "p")],
            *[f"v_{n}_pp{i}{t}" for n in ("2", "5", "6") for i in ("s", "p") for t in ("", "_tb")],
            *[f"v_{n}_dd{i}" for n in ("2", "5", "6") for i in ("s", "p", "d")]
        ]

        sk_simple_params_names = [
            *[rf"$V^{{{n}}}_{{pd{i}}}$" for n in ("1", "3", "4") for i in (r"\sigma", r"\pi")],
            *[rf"$V^{{{n}}}_{{pp{i}{t}}}$" for n in ("2", "5", "6") for i in (r"\sigma", r"\pi") for t in ("", ",tb")],
            *[rf"$V^{{{n}}}_{{dd{i}}}$" for n in ("2", "5", "6") for i in (r"\sigmas", r"\pi", r"\delta")]
        ]
        self._sk_simple_params_dict: dict = dict(zip(
            sk_simple_params,
            [FloatParameter(name=p_name) for p_name in sk_simple_params_names]
        ))
        self._sk_params_changable = [
            "theta",
            *[f"delta_{r}" for r in ("p", "z", "0", "1", "2")],
            *[f"v_0_pp{r}" for r in ("s", "p")]
        ]
        if input_dict is not None:
            self.from_dict(input_dict)

    def _recalculate_params(self):
        if self._recalculate_params_bool:
            # h1
            self._sk_params_dict["v_1_e_pds"].param = self["v_1_pds"]
            self._sk_params_dict["v_1_o_pds"].param = self["v_1_pds"]
            self._sk_params_dict["v_1_e_pdp"].param = self["v_1_pdp"]
            self._sk_params_dict["v_1_o_pdp"].param = self["v_1_pdp"]
            # h2
            self._sk_params_dict["v_2_e_dds"].param = self["v_2_dds"]
            self._sk_params_dict["v_2_e_ddp"].param = self["v_2_ddp"]
            self._sk_params_dict["v_2_o_ddp"].param = self["v_2_ddp"]
            self._sk_params_dict["v_2_e_ddd"].param = self["v_2_ddd"]
            self._sk_params_dict["v_2_o_ddd"].param = self["v_2_ddd"]
            self._sk_params_dict["v_2_e_pps"].param = self["v_2_pps"]
            self._sk_params_dict["v_2_o_pps"].param = self["v_2_pps"]
            self._sk_params_dict["v_2_e_ppp"].param = self["v_2_ppp"]
            self._sk_params_dict["v_2_o_ppp"].param = self["v_2_ppp"]
            self._sk_params_dict["v_2_e_pps_tb"].param = self["v_2_pps_tb"]
            self._sk_params_dict["v_2_o_pps_tb"].param = self["v_2_pps_tb"]
            self._sk_params_dict["v_2_e_ppp_tb"].param = self["v_2_ppp_tb"]
            self._sk_params_dict["v_2_o_ppp_tb"].param = self["v_2_ppp_tb"]
            # h3
            self._sk_params_dict["v_3_e_pds"].param = self["v_3_pds"]
            self._sk_params_dict["v_3_o_pds"].param = self["v_3_pds"]
            self._sk_params_dict["v_3_e_pdp"].param = self["v_3_pdp"]
            self._sk_params_dict["v_3_o_pdp"].param = self["v_3_pdp"]
            # h4
            self._sk_params_dict["v_4_e_pds"].param = self["v_4_pds"]
            self._sk_params_dict["v_4_o_pds"].param = self["v_4_pds"]
            self._sk_params_dict["v_4_e_pdp"].param = self["v_4_pdp"]
            self._sk_params_dict["v_4_o_pdp"].param = self["v_4_pdp"]
            # h5
            self._sk_params_dict["v_5_e_dds"].param = self["v_5_dds"]
            self._sk_params_dict["v_5_e_ddp"].param = self["v_5_ddp"]
            self._sk_params_dict["v_5_o_ddp"].param = self["v_5_ddp"]
            self._sk_params_dict["v_5_e_ddd"].param = self["v_5_ddd"]
            self._sk_params_dict["v_5_o_ddd"].param = self["v_5_ddd"]
            self._sk_params_dict["v_5_e_pps"].param = self["v_5_pps"]
            self._sk_params_dict["v_5_o_pps"].param = self["v_5_pps"]
            self._sk_params_dict["v_5_e_ppp"].param = self["v_5_ppp"]
            self._sk_params_dict["v_5_o_ppp"].param = self["v_5_ppp"]
            self._sk_params_dict["v_5_e_pps_tb"].param = self["v_5_pps_tb"]
            self._sk_params_dict["v_5_o_pps_tb"].param = self["v_5_pps_tb"]
            self._sk_params_dict["v_5_e_ppp_tb"].param = self["v_5_ppp_tb"]
            self._sk_params_dict["v_5_o_ppp_tb"].param = self["v_5_ppp_tb"]
            # h6
            self._sk_params_dict["v_6_e_dds"].param = self["v_6_dds"]
            self._sk_params_dict["v_6_e_ddp"].param = self["v_6_ddp"]
            self._sk_params_dict["v_6_o_ddp"].param = self["v_6_ddp"]
            self._sk_params_dict["v_6_e_ddd"].param = self["v_6_ddd"]
            self._sk_params_dict["v_6_o_ddd"].param = self["v_6_ddd"]
            self._sk_params_dict["v_6_e_pps"].param = self["v_6_pps"]
            self._sk_params_dict["v_6_o_pps"].param = self["v_6_pps"]
            self._sk_params_dict["v_6_e_ppp"].param = self["v_6_ppp"]
            self._sk_params_dict["v_6_o_ppp"].param = self["v_6_ppp"]
            self._sk_params_dict["v_6_e_pps_tb"].param = self["v_6_pps_tb"]
            self._sk_params_dict["v_6_o_pps_tb"].param = self["v_6_pps_tb"]
            self._sk_params_dict["v_6_e_ppp_tb"].param = self["v_6_ppp_tb"]
            self._sk_params_dict["v_6_o_ppp_tb"].param = self["v_6_ppp_tb"]
            super()._recalculate_params()

    @property
    def _allowed_params(self):
        return [*self._general_params_dict, *self._energy_params_dict, *self._sk_params_dict, *self._sk_simple_params_dict]

    def __setitem__(self, key, value):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                self._general_params_dict[key].param = value
            elif key in [*self._sk_simple_params_dict]:
                self._sk_simple_params_dict[key].param = value
                self._recalculate_params()
            elif key in [*self._sk_params_dict]:
                if not (key in self._sk_params_changable):
                    warnings.warn("The variable {0} is set by SK functions and is read-only, don't change it.".format(key))
                self._sk_params_dict[key].param = value
                self._recalculate_params()
            elif key in [*self._energy_params_dict]:
                warnings.warn("The variable {0} is set by SK functions and is read-only, don't change it.".format(key))
                self._energy_params_dict[key].param = value
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def get_param(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].param
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].param
            elif key in [*self._sk_params_dict]:
                return self._sk_params_dict[key].param
            elif key in [*self._sk_simple_params_dict]:
                return self._sk_simple_params_dict[key].param
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

    def get_name(self, key):
        if self._check_key(key):
            if key in [*self._general_params_dict]:
                return self._general_params_dict[key].name
            elif key in [*self._energy_params_dict]:
                return self._energy_params_dict[key].name
            elif key in [*self._sk_params_dict]:
                return self._sk_params_dict[key].name
            elif key in [*self._sk_simple_params_dict]:
                return self._sk_simple_params_dict[key].name
            else:
                warnings.warn("This should not happen, {0} should be a valid key.".format(key))

