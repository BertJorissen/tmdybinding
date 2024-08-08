from typing import Optional, Dict, Union, List
from .symmetry_group import ParametersList, FloatParameter

class StrainParametersList(ParametersList):
    """List of strain parameters."""

    def __init__(self, input_dict: Optional[Dict[str, Union[float, List[float]]]] = None, max_order: int = 1):
        super().__init__(None)
        self.max_order = max_order
        strain_params = []
        strain_params_names = []

        # biaxial and uniaxial strain parameters
        for strain_type in ("b", "u", "s"):
            for order_idx in range(1, self.max_order+1):
                strain_part = f"{strain_type}_{str(order_idx)}"
                strain_part_name = rf"{strain_type},{str(order_idx)}"
                if not strain_type == "s":
                    for param_key, param_float in self._energy_params_dict.items():
                        strain_params.append(rf"{param_key}_{strain_part}")
                        check_for = r"$\epsilon_" if "eps_" in param_key else r"$u_"
                        split_idx, split_len = len(check_for), 1
                        param_name = param_float.name
                        pre_name = param_name[:split_idx]
                        idx_name = param_name[split_idx:split_idx+split_len]
                        pst_name = param_name[split_idx+split_len:]
                        assert pre_name == check_for, f"Expected {check_for}, got {pre_name}"
                        strain_params_names.append(rf"{pre_name}{{{idx_name},{strain_part_name}}}{pst_name}")
                else:
                    params_shear = [
                        *[f"u_{u}_{i}_m_e_{strain_part}" for u in ("1", "3", "4") for i in range(5, 9)],
                        *[f"u_{u}_{i}_m_o_{strain_part}" for u in ("1", "3", "4") for i in range(3, 6)],
                        *[f"u_{u}_{i}_{m}_e_{strain_part}" for u in ("2", "6") for i in ("1", "2", "4") for m in ("m", "x")],
                        *[f"u_{u}_{i}_x_o_{strain_part}" for u in ("2", "6") for i in ("1", "2", "4")],
                        *[f"u_{u}_{i}_m_o_{strain_part}" for u in ("2", "6") for i in range(1)],
                        *[f"u_{u}_{i}_m_e_{strain_part}" for u in "5" for i in ("2", "4", "7", "8")],
                        *[f"u_{u}_{i}_x_{r}_{strain_part}" for u in "5" for i in ("1", "4", "7", "8") for r in ("e", "o")],
                        *[f"u_{u}_{i}_m_o_{strain_part}" for u in "5" for i in ("1", "3")],
                    ]
                    params_names_shear = [
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},e}}$" for u in ("1", "3", "4") for i in range(5, 9)],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},o}}$" for u in ("1", "3", "4") for i in range(3, 6)],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},{m}e}}$" for u in ("2", "6") for i in ("1", "2", "4") for m in ("M", "X")],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},Xo}}$" for u in ("2", "6") for i in ("1", "2", "4")],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},Mo}}$" for u in ("2", "6") for i in range(1)],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},Me}}$" for u in "5" for i in ("2", "4", "7", "8")],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},X{r}}}$" for u in "5" for i in ("2", "4", "7", "8") for r in ("e", "o")],
                        *[rf"$u_{{{u},{strain_part_name}}}^{{{i},Mo}}$" for u in "5" for i in ("1", "3")],
                    ]
                    strain_params.extend(params_shear)
                    strain_params_names.extend(params_names_shear)


        # shear strain parameters, only hoppings as onsite are given by the uniaxial parameters
        self._strain_params_dict: dict = dict(zip(
            strain_params,
            [FloatParameter(name=p_name) for p_name in strain_params_names]
        ))


        if input_dict is not None:
            self.from_dict(input_dict)

    @property
    def _unique_params_dict(self) -> List[dict]:
        return [self._general_params_dict, self._energy_params_dict, self._strain_params_dict]

