import numpy as np
from typing import Optional, Dict, Union, List
from .symmetry_group import ParametersList, FloatParameter


class SKParametersList(ParametersList):
    """Class to store the parameters for SK e/o model"""

    def __init__(self, input_dict: Optional[Dict[str, Union[str, Optional[float]]]] = None):
        """Initialize the SK parameters

        Parameters:
            input_dict (dict): The dictionary containing the SK parameters. The keys are the parameter names and the
                values are the parameter values. The keys are:
                `a`, `lamb_x`, `lamb_m`, `material`, `theta`,
                `delta_p`, `delta_z`, `delta_0`, `delta_1`, `delta_2`,
                `v_0_pps`, `v_0_ppp`,
                `v_1_e_pds`, `v_1_e_pdp`, `v_1_o_pds`, `v_1_o_pdp`,
                `v_2_e_pps`, `v_2_e_ppp`, `v_2_e_pps_tb`, `v_2_e_ppp_tb`,
                `v_2_o_pps`, `v_2_o_ppp`, `v_2_o_pps_tb`, `v_2_o_ppp_tb`,
                `v_2_e_dds`, `v_2_e_ddp`, `v_2_e_ddd`,
                `v_2_o_ddp`, `v_2_o_ddd`,
                `v_3_e_pds`, `v_3_e_pdp`, `v_3_o_pds`, `v_3_o_pdp`,
                `v_4_e_pds`, `v_4_e_pdp`, `v_4_o_pds`, `v_4_o_pdp`,
                `v_5_e_pps`, `v_5_e_ppp`, `v_5_e_pps_tb`, `v_5_e_ppp_tb`,
                `v_5_o_pps`, `v_5_o_ppp`, `v_5_o_pps_tb`, `v_5_o_ppp_tb`,
                `v_5_e_dds`, `v_5_e_ddp`, `v_5_e_ddd`,
                `v_5_o_ddp`, `v_5_o_ddd`,
                `v_6_e_pps`, `v_6_e_ppp`, `v_6_e_pps_tb`, `v_6_e_ppp_tb`,
                `v_6_o_pps`, `v_6_o_ppp`, `v_6_o_pps_tb`, `v_6_o_ppp_tb`,
                `v_6_e_dds`, `v_6_e_ddp`, `v_6_e_ddd`,
                `v_6_o_ddp` and `v_6_o_ddd`"""
        super().__init__(None)
        self._recalculate_params_bool: bool = True
        self._tan_theta: Optional[float] = None

        onsite_sk_params = [
            "theta",
            *[f"delta_{r}" for r in ("p", "z", "0", "1", "2")],
            *[f"v_0_pp{r}" for r in ("s", "p")]
        ]
        onsite_sk_params_names = [
            r"$\theta$",
            *[rf"$\delta_{r}$" for r in ("p", "z", "0", "1", "2")],
            *[rf"$V^0_{{pp{r}}}$" for r in (r"\sigma", r"\pi")],
        ]

        self._onsite_sk_params_dict: dict = dict(zip(
            onsite_sk_params,
            [FloatParameter(name=p_name) for p_name in onsite_sk_params_names]
        ))

        sk_params = [
            *[f"v_{n}_{r}_pd{i}" for n in ("1", "3", "4") for r in ("e", "o") for i in ("s", "p")],
            *[f"v_{n}_{r}_pp{i}{t}" for n in ("2", "5", "6")
              for r in ("e", "o") for i in ("s", "p") for t in ("", "_tb")],
            *[f"v_{n}_e_dds" for n in ("2", "5", "6")],
            *[f"v_{n}_{r}_dd{i}" for n in ("2", "5", "6") for r in ("e", "o") for i in ("p", "d")]
        ]

        sk_params_names = [
            *[rf"$V^{{{n}_{r}}}_{{pd{i}}}$" for n in ("1", "3", "4") for r in ("e", "o") for i in (r"\sigma", r"\pi")],
            *[rf"$V^{{{n}_{r}}}_{{pp{i}{t}}}$" for n in ("2", "5", "6")
              for r in ("e", "o") for i in (r"\sigma", r"\pi") for t in ("", ",tb")],
            *[rf"$V^{{{n}_e}}_{{dd\sigmas}}$" for n in ("2", "5", "6")],
            *[rf"$V^{{{n}_{r}}}_{{dd{i}}}$" for n in ("2", "5", "6") for r in ("e", "o") for i in (r"\pi", r"\delta")]
        ]
        self._sk_params_dict: dict = dict(zip(
            sk_params,
            [FloatParameter(name=p_name) for p_name in sk_params_names]
        ))

        if input_dict is not None:
            self.from_dict(input_dict)

    @property
    def _unique_params_dict(self) -> List[dict]:
        return [self._general_params_dict, self._onsite_sk_params_dict, self._sk_params_dict]

    @property
    def _protected_params_dict(self) -> Optional[List[dict]]:
        return [self._energy_params_dict]

    def __setitem__(self, key, value):
        """Set the value of a parameter and recalculate the parameters."""
        super().__setitem__(key, value)
        self._recalculate_params()

    def from_dict(self, input_dict: Dict[str, Union[float, str]]):
        """Set the parameters from a dictionary and recalculate the parameters.

        Parameters:
            input_dict (dict): The dictionary containing the SK parameters. The keys are the parameter names and the
                values are the parameter values. The keys are:
                `a`, `lamb_x`, `lamb_m`, `material`, `theta`,
                `delta_p`, `delta_z`, `delta_0`, `delta_1`, `delta_2`,
                `v_0_pps`, `v_0_ppp`,
                `v_1_e_pds`, `v_1_e_pdp`, `v_1_o_pds`, `v_1_o_pdp`,
                `v_2_e_pps`, `v_2_e_ppp`, `v_2_e_pps_tb`, `v_2_e_ppp_tb`,
                `v_2_o_pps`, `v_2_o_ppp`, `v_2_o_pps_tb`, `v_2_o_ppp_tb`,
                `v_2_e_dds`, `v_2_e_ddp`, `v_2_e_ddd`,
                `v_2_o_ddp`, `v_2_o_ddd`,
                `v_3_e_pds`, `v_3_e_pdp`, `v_3_o_pds`, `v_3_o_pdp`,
                `v_4_e_pds`, `v_4_e_pdp`, `v_4_o_pds`, `v_4_o_pdp`,
                `v_5_e_pps`, `v_5_e_ppp`, `v_5_e_pps_tb`, `v_5_e_ppp_tb`,
                `v_5_o_pps`, `v_5_o_ppp`, `v_5_o_pps_tb`, `v_5_o_ppp_tb`,
                `v_5_e_dds`, `v_5_e_ddp`, `v_5_e_ddd`,
                `v_5_o_ddp`, `v_5_o_ddd`,
                `v_6_e_pps`, `v_6_e_ppp`, `v_6_e_pps_tb`, `v_6_e_ppp_tb`,
                `v_6_o_pps`, `v_6_o_ppp`, `v_6_o_pps_tb`, `v_6_o_ppp_tb`,
                `v_6_e_dds`, `v_6_e_ddp`, `v_6_e_ddd`,
                `v_6_o_ddp` and `v_6_o_ddd`"""
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
        u_0_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_1: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_4: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = -vppstb if vstb_bool else None
        u_5_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
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
        u_0_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = -vppstb if vstb_bool else None
        u_5_2: Optional[float] = -r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        u_6: Optional[float] = 2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        return [u_0, u_2, u_3, u_5, u_6]

    def _h_xx_x_o(self, r: float, vpps: Optional[float], vppp: Optional[float],
                  vppstb: Optional[float] = None, vppptb: Optional[float] = None) -> List[Optional[float]]:
        vpps_bool = vpps is not None
        vppp_bool = vppp is not None
        vptb_bool = vppstb is not None
        vstb_bool = vppptb is not None
        u_0_0: Optional[float] = vpps if vpps_bool else None
        u_0_1: Optional[float] = -vppptb if vptb_bool else None
        u_0_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_1: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_2: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = -vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_4: Optional[float] = 0 if (vpps_bool or vppp_bool or vstb_bool or vptb_bool) else None
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = vppstb if vstb_bool else None
        u_5_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
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
        u_0_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_0: Optional[float] = self._comb_nonefloat([u_0_0, u_0_1, u_0_2])
        u_2: Optional[float] = 2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_3_0: Optional[float] = vppp if vppp_bool else None
        u_3_1: Optional[float] = -vppptb if vptb_bool else None
        u_3: Optional[float] = self._comb_nonefloat([u_3_0, u_3_1])
        u_5_0: Optional[float] = vppp if vppp_bool else None
        u_5_1: Optional[float] = vppstb if vstb_bool else None
        u_5_2: Optional[float] = r ** 2 * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        u_5: Optional[float] = self._comb_nonefloat([u_5_0, u_5_1, u_5_2])
        u_6: Optional[float] = -2 * self._tan_theta * r * (vppptb - vppstb) / (r ** 2 + 4 * self._tan_theta ** 2)\
            if vstb_bool and vptb_bool else None
        return [u_0, u_2, u_3, u_5, u_6]

    def _h_mx_e(self, r: float, vpds: Optional[float], vpdp: Optional[float]) -> List[Optional[float]]:
        if (vpds is not None) and (vpdp is not None):
            u_0 = np.sqrt(2) * r * vpdp / (np.sqrt(r ** 2 + self._tan_theta ** 2))
            u_1 = -r * (
                    2 * np.sqrt(3) * self._tan_theta ** 2 * vpdp + (r ** 2 - 2 * self._tan_theta ** 2) * vpds
            ) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_2 = -r * (
                    2 * self._tan_theta ** 2 * vpdp + np.sqrt(3) * r ** 2 * vpds
            ) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_3 = self._tan_theta * (
                    r ** 2 * 2 * np.sqrt(3) * vpdp + (2 * self._tan_theta ** 2 - r ** 2) * vpds
            ) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            u_4 = r ** 2 * self._tan_theta * (
                    2 * vpdp - np.sqrt(3) * vpds
            ) / (np.sqrt(2) * ((r ** 2 + self._tan_theta ** 2) ** (3 / 2)))
            return [u_0, u_1, u_2, u_3, u_4]
        else:
            return [None, None, None, None, None]

    def _h_mx_o(self, r: float, vpds: Optional[float], vpdp: Optional[float]) -> List[Optional[float]]:
        if (vpds is not None) and (vpdp is not None):
            u_0 = np.sqrt(2) * self._tan_theta * vpdp / np.sqrt(r ** 2 + self._tan_theta ** 2)
            u_1 = np.sqrt(2) * self._tan_theta * (
                    (self._tan_theta ** 2 - r ** 2) * vpdp + r ** 2 * np.sqrt(3) * vpds
            ) / ((r ** 2 + self._tan_theta ** 2) ** (3 / 2))
            u_2 = np.sqrt(2) * r * (
                    (r ** 2 - self._tan_theta ** 2) * vpdp + self._tan_theta ** 2 * np.sqrt(3) * vpds
            ) / ((r ** 2 + self._tan_theta ** 2) ** (3 / 2))
            return [u_0, u_1, u_2]
        else:
            return [None, None, None]


class SKSimpleParametersList(SKParametersList):
    """Class to store the simplified SK paramters (no even/odd)."""
    def __init__(self, input_dict: Optional[Dict[str, Union[str, Optional[float]]]] = None):
        """Initializes the SKSimpleParametersList object.

        Parameters:
            input_dict (dict): Dictionary with the SK paramters. The keus are the parameter names and the values are the
                parameter values. The parameter names are:
                `a`, `material`, `lamb_x`, `lamb_m`, `theta`,
                `delta_p`, `delta_z`, `delta_0`, `delta_1`, `delta_2`,
                `v_0_pps`, `v_0_ppp`,
                `v_1_pds`, `v_1_pdp`,
                `v_2_dds`, `v_2_ddp`, `v_2_ddd`, `v_2_pps`, `v_2_ppp`, `v_2_pps_tb`, `v_2_ppp_tb`,
                `v_3_pds`, `v_3_pdp`,
                `v_4_pds`, `v_4_pdp`,
                `v_5_dds`, `v_5_ddp`, `v_5_ddd`, `v_5_pps`, `v_5_ppp`, `v_5_pps_tb`, `v_5_ppp_tb`,
                `v_6_dds`, `v_6_ddp`, `v_6_ddd`, `v_6_pps`, `v_6_ppp`, `v_6_pps_tb` and `v_6_ppp_tb`."""
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
    def _unique_params_dict(self) -> List[dict]:
        return [self._general_params_dict, self._onsite_sk_params_dict, self._sk_simple_params_dict]

    @property
    def _protected_params_dict(self) -> Optional[List[dict]]:
        return [self._sk_params_dict, self._energy_params_dict]


_params_sk_dias_mos2_fitted = SKParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.31882931,     0.0806,         0.0536,         "MoS2",         0.703,
        -6.32462197e+00, -4.77732844e+00, -5.90160369e+00, -9.523165605, -12.5025582,
        5.630393, -2.447044695,
        4.00851699e+00, -1.59507645e+00,        2.15787943e+00, -1.06942019e+00,
        -6.58197057e-01,  5.72293847e-01,  2.61828424e-01,  9.92502370e-01, -1.66848188e+00,
        4.27331648e-02,  1.42165455e-02, 5.74382013e-01, -1.53696505e-01,
        3.22763559e-02,  5.67853736e-02, -3.79415525e-02,  1.60731513e-01,  1.50517193e-03,
        7.10141101e-04, 6.55934811e-02, -3.36741023e-02, -2.23077491e-02
    ]
)))

_params_sk_rostami_fitted = SKSimpleParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.31882931,     0.075,          0.052,          "MoS2",         0.703,
        -6.2539434,     0.000,          -5.64168266,    -10.24784771,   -16.61743716,
        3.89201485,     -1.36521723,    -0.73264736,    0.65489363,     0.26205299,
        0.000,          0.000,          0.73420676,     -1.44859604,
    ]
)))

_params_sk_cappelluti_fitted = SKSimpleParametersList(dict(zip(
    # fitted https://doi.org/10.21468/SciPostPhysCore.7.1.004
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.31882931,     0.075,          0.052,          "MoS2",         0.703,
        -6.2539434,     -7.6674323,     -5.64168266,    -10.24784771,   -9.10842771,
        3.89201485,     -1.36521723,    -0.73264736,    0.65489363,     0.26205299,
        7.50900945,     -0.62850624,    0.73420676,     -1.44859604,
    ]
)))

_params_sk_cappelluti = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.88.075409
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.512,         None,           -3.025,         -1.276,         -8.236,
        -2.619,         -1.396,         -0.933,         -0.478,         -0.442,
        0.696,          0.278,          0.696,           0.278,
    ]
)))

_params_sk_roldan_mos2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/2053-1583/1/3/034003
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.512,         0.419,          -3.025,         -1.276,         -8.236,
        -2.619,         -1.396,         -0.933,         -0.478,         -0.442,
        0.696,          0.278,          0.696,          0.278
    ]
)))

_params_sk_roldan_ws2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/2053-1583/1/3/034003
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3153,         0.215,          0.057,          "WS2",          0.712,
        -1.550,         0.851,          -3.090,         -1.176,         -7.836,
        -0.619,         -1.396,         -0.983,         -0.478,         -0.442,
        0.696,          0.278,          0.696,          0.278
    ]
)))

_params_sk_ridolfi = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        0.201,          -1.563,         -0.352,         -54.839,        -39.275,
        -9.88,          4.196,          -1.153,         0.612,          0.086,
        12.734,         -2.175,         12.734,         -2.175
    ]
)))


_params_sk_ridolfi_vb = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        0.191,          -1.599,         0.081,          -48.934,        -37.981,
        -8.963,         4.115,          -1.154,         0.964,          0.117,
        10.707,         -4.084,         10.707,         -4.084
    ]
)))

_params_sk_ridolfi_minimal = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        -11.683,        -208.435,       -75.952,        -23.761,        -35.968,
        -56.738,        1.318,          -2.652,         1.750,          1.482,
        0.000,          0.000,          0.000,          0.000
    ]
)))

_params_sk_venkateswarlu = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.102.081103
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.318,          0.086,          0.052,          "MoS2",         0.704,
        0.1356,         -0.4204,        0.0149,         -38.71,         -29.45,
        -7.193,         3.267,          -0.9035,         0.7027,        0.0897,
        8.079,          -2.678,         7.336,          -2.432
    ]
)))

_params_sk_silva_guillen_mos2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.3390/app6100284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.086,          0.052,          "MoS2",         0.716,
        -1.094,         -0.050,         -1.511,         -3.559,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))

_params_sk_silva_guillen_mose2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.3390/app6100284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.288,          0.089,          0.256,          "MoSe2",        0.720,
        -1.144,         -0.250,         -1.488,         -4.931,         -7.503,
        3.728,          -1.222,         -0.823,         0.215,          0.192,
        1.256,          -0.205,         1.256,          -0.205
    ]
)))

_params_sk_silva_guillen_ws2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.3390/app6100284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.153,          0.271,          0.057,          "WS2",          0.712,
        -1.155,         -0.650,         -2.279,         -3.864,         -7.327,
        7.911,          -1.220,         -1.328,         0.121,          0.422,
        1.178,          -0.273,         1.178,          -0.273
    ]
)))

_params_sk_silva_guillen_wse2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.3390/app6100284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.260,          0.251,          0.439,          "WSe2",         0.722,
        -0.935,         -1.250,         -2.321,         -5.629,         -6.759,
        5.803,          -1.081,         -1.129,         0.094,          0.317,
        1.530,          -0.123,         1.530,          -0.123
    ]
)))

_params_sk_pearce = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.94.155416
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.319,          0.075,          0.052,          "MoS2",         0.687,
        2.12,           -0.46,          -1.41,          -5.38,          -3.69,
        -2.83,          0.67,           -0.24,          -0.62,          0.45,
        -0.42,          -1.32,          -0.42,          -1.32
    ]
)))

_params_sk_bieniek = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.97.085153
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3193,         0.075,          0.052,          "MoS2",         0.71,
        -0.03,          -0.03,          -0.03,          -3.36,          -4.78,
        -3.39,          1.1,            -1.1,           0.76,           0.27,
        1.19,           -0.83,          1.19,           -0.83
    ]
)))

_params_sk_abdi_mos2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3166,         0.086,          0.052,          "MoS2",         0.710,
        -1.094,         -0.050,         -1.511,         -3.559,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))

_params_sk_abdi_wse2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3260,         0.251,          0.439,          "WSe2",         0.722,
        -0.935,         -1.250,         -2.321,         -5.629,         -6.759,
        5.803,          -1.081,         -1.129,         0.094,          0.317,
        1.530,          -0.123,         1.530,          -0.205
    ]
)))

_params_sk_abdi_mose2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3298,         0.089,          0.256,          "MoSe2",        0.721,
        -1.144,         -0.250,         -1.488,         -4.931,         -7.503,
        3.728,          -1.222,         -0.823,         0.215,          0.192,
        1.256,          -0.205,         1.256,          -0.205,
    ]
)))

_params_sk_rostami = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.92.195402
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.094,         None,           -1.512,         -3.560,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))

_params_sk_dias_mos2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3166,         0.0806,         0.0536,         "MoS2",         0.710,
        -0.4939,        0.5624,         -0.2473,        -3.14115,       -4.0595,
        4.1989,         -1.89235,
        4.2398,         -1.2413,        2.2251,         -0.7614,
        -0.6717,        0.5706,         0.2729,         -0.0914,        -0.4619,
        0.0150,         0.0497,         0.8131,         -0.2763,
        0.0314,         0.0961,         -0.0305,        0.3723,         0.0014,
        0.0051,         0.0184,         -0.0395,        0.0092
    ]
)))

_params_sk_dias_mose2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3288,         0.0806,         0.0820,         "MoSe2",        0.710,
        -0.1276,        0.3046,         -0.2724,        -3.8352,        -4.30195,
        4.30095,        -2.8093,
        3.4524,         -1.4295,        2.0197,         -0.6811,
        -0.6674,        0.5573,         0.0970,         1.2630,         -0.4857,
        0.01637,        0.0965,         0.9449,         -0.3039,
        0.0776,         0.0573,         -0.04778,       0.2372,         0.0249,
        0.0140,         0.0354,         -0.0293,        -0.0094
    ]
)))

_params_sk_dias_mote2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3519,         0.0806,         0.1020,         "MoTe2",        0.710,
        -0.6630,        0.0491,         -0.2852,        -0.9084,        -1.8434,
        2.6799,         0.0678,
        2.2362,         -0.6279,        1.8294,         -0.5048,
        -0.4795,        -0.0934,        0.1656,         0.8198,         -0.2483,
        0.3267,         0.3033,         0.8459,         -0.4143,
        -0.1493,        -0.0627,        0.0360,         0.1169,         0.2683,
        -0.0617,        0.1002,         0.0114,         -0.0092
    ]
)))

_params_sk_dias_ws2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.31532,        0.2754,         0.0536,         "WS2",          0.710,
        -0.3609,        0.8877,         -0.7364,        -3.52825,       -4.5926,
        4.415,          -1.97685,
        5.2769,         -1.2119,        2.4044,         -0.8115,
        -0.8942,        0.7347,         0.3417,         -0.3943,        -0.4069,
        -0.0142,        0.0036,         0.8415,         -0.2661,
        0.0508,         0.1278,         -0.0091,        0.1415,         0.0261,
        -0.0135,        -0.0191,        -0.0169,        0.0262
    ]
)))

_params_sk_dias_wse2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3282,         0.2754,         0.0820,         "WSe2",         0.710,
        -0.5558,        0.6233,         -1.9340,        -2.20345,       -3.6177,
        3.1056,         -0.99385,
        5.1750,         -0.9139,        2.1733,         -0.7688,
        -0.8697,        0.6206,         0.3743,         0.1311,         -0.2475,
        -0.0469,        0.0923,         0.9703,         -0.2920,
        0.0443,         0.0912,         -0.0447,        0.1197,         0.1075,
        0.0096,         0.0140,         -0.0451,        0.0113
    ]
)))

_params_sk_peng6_mos2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3160,         0.074,          0.052,          "MoS2",         0.716,
        -1.094,         0.000,          -1.511,         -3.559,         -6.886,
        -3.679,         1.199,          -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))


_params_sk_peng6_ws2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3153,         0.210,          0.159,          "WS2",          0.712,
        -1.090,         0.000,          -1.525,         -3.599,         -7.598,
        -3.785,         1.275,          -0.925,         0.261,          0.220,
        1.250,          -0.476,         1.250,          -0.476
    ]
)))

_params_sk_peng11_mos2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3160,         0.1050,         0.0536,         "MoS2",         0.716,
        -0.5048,        0.5624,         -0.2399,        -3.1612,        -4.5063,
        4.7433,         -1.9457,
        4.9182,         -1.4399,        2.5811,         -0.8832,
        -0.8060,        0.6847,         0.3275,         -0.1024,        -0.5173,
        0.0180,         0.0596,         0.9107,         -0.3095,
        0.0377,         0.1153,         -0.0366,        0.4170,         0.0016,
        0.0061,         0.0221,         -0.0442,        0.0103
    ]
)))

_params_sk_peng11_ws2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3153,         0.3188,         0.0536,         "WS2",          0.717,
        -0.3429,        0.9765,         -0.8690,        -2.8784,        -5.1765,
        5.0999,         -2.2677,
        6.1212,         -1.4058,        2.7891,         -0.9413,
        -1.0730,        0.8816,         0.41004,        -0.4416,        -0.4557,
        -0.0170,        0.00432,        0.9425,         -0.2980,
        0.0609,         0.1534,         -0.0109,        0.1585,         0.0292,
        -0.0162,        -0.0229,        -0.0189,        0.0293
    ]
)))

_params_sk_peng11_mose2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3288,         0.1050,         0.1680,         "MoSe2",        0.714,
        -0.2080,        1.1248,         -0.0840,        -3.2502,        -3.3814,
        3.6329,         -1.6049,
        4.3280,         -1.2671,        2.2714,         -0.7772,
        -0.9269,        0.7874,         0.2784,         -0.0932,        -0.5639,
        0.0207,         0.0507,         0.9926,         -0.3373,
        0.0320,         0.1326,         -0.0421,        0.3795,         0.0014,
        0.0070,         0.0188,         -0.0403,        0.0094
    ]
)))

_params_sk_peng11_wse2 = SKParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.109.245412
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3282,         0.3188,         0.1680,         "WSe2",         0.714,
        -0.3799,        0.5219,         -0.9623,        -3.6513,        -4.2498,
        3.3328,         -2.2857,
        5.4356,         -1.2484,        2.4767,         -0.8359,
        -1.2233,        1.0051,         0.3526,         -0.4045,        -0.4940,
        -0.0147,        0.0037,         1.0216,         -0.3231,
        0.0524,         0.1748,         -0.0124,        0.1718,         0.0268,
        -0.0139,        -0.0261,        -0.0173,        0.0269
    ]
)))

rostami = {
    "MoS2": _params_sk_rostami,
    "fitted": {"MoS2": _params_sk_rostami_fitted}
}

dias = {
    "MoS2": _params_sk_dias_mos2,
    "MoSe2": _params_sk_dias_mose2,
    "MoTe2": _params_sk_dias_mote2,
    "WS2": _params_sk_dias_ws2,
    "WSe2": _params_sk_dias_wse2,
    "fitted": {"MoS2": _params_sk_dias_mos2_fitted}
}

cappelluti = {
    "MoS2": _params_sk_cappelluti,
    "fitted": {"MoS2": _params_sk_cappelluti_fitted}
}

roldan = {
    "MoS2": _params_sk_roldan_mos2,
    "WS2": _params_sk_roldan_ws2
}

ridolfi = {
    "MoS2": {
        "normal": _params_sk_ridolfi,
        "vb": _params_sk_ridolfi_vb,
        "minimal": _params_sk_ridolfi_minimal
    }
}

venkateswarlu = {
    "MoS2": _params_sk_venkateswarlu
}

silva_guillen = {
    "MoS2": _params_sk_silva_guillen_mos2,
    "MoSe2": _params_sk_silva_guillen_mose2,
    "WS2": _params_sk_silva_guillen_ws2,
    "WSe2": _params_sk_silva_guillen_wse2
}

pearce = {
    "MoS2": _params_sk_pearce
}

bieniek = {
    "MoS2": _params_sk_bieniek
}

abdi = {
    "MoS2": _params_sk_abdi_mos2,
    "MoSe2": _params_sk_abdi_mose2,
    "WS2": _params_sk_abdi_wse2
}

peng6 = {
    "MoS2": _params_sk_peng6_mos2,
    "WS2": _params_sk_peng6_ws2
}

peng11 = {
    "MoS2": _params_sk_peng11_mos2,
    "MoSe2": _params_sk_peng11_mose2,
    "WS2": _params_sk_peng11_ws2,
    "WSe2": _params_sk_peng11_wse2
}
