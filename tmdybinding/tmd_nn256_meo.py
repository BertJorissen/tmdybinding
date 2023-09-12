"""Tight-binding models for group 1 transition metal dichalcogenides (tmd), 5 band."""
from .tmd_abstract_lattice import AbstractLattice
from .tmd_matrices import TmdMatrices


class TmdNN256Meo(AbstractLattice):
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
