#
# Diffusion-limited kinetics
#

import pybamm
from ..base_interface import BaseInterface


class DiffusionLimited(BaseInterface):
    """
    Submodel for diffusion-limited kinetics

    Parameters
    ----------
    param :
        model parameters
    domain : str
        The domain to implement the model, either: 'Negative' or 'Positive'.
    reaction : str
        The name of the reaction being implemented
    options: dict
        A dictionary of options to be passed to the model. See
        :class:`pybamm.BaseBatteryModel`
    order : str
        The order of the model ("leading" or "full")

    **Extends:** :class:`pybamm.interface.BaseInterface`
    """

    def __init__(self, param, domain, reaction, options, order):
        super().__init__(param, domain, reaction, options)
        self.order = order

    def get_coupled_variables(self, variables):
        domain, Domain = self.domain_Domain
        reaction_name = self.reaction_name

        delta_phi_s = variables[f"{Domain} electrode surface potential difference"]
        # If delta_phi_s was broadcast, take only the orphan
        if isinstance(delta_phi_s, pybamm.Broadcast):
            delta_phi_s = delta_phi_s.orphans[0]

        # Get exchange-current density
        j0 = self._get_exchange_current_density(variables)
        # Get open-circuit potential variables and reaction overpotential
        ocp = variables[f"{Domain} electrode {reaction_name}open circuit potential"]
        eta_r = delta_phi_s - ocp

        # Get interfacial current densities
        j = self._get_diffusion_limited_current_density(variables)
        j_tot_av, a_j_tot_av = self._get_average_total_interfacial_current_density(
            variables
        )

        variables.update(self._get_standard_interfacial_current_variables(j))
        variables.update(
            self._get_standard_total_interfacial_current_variables(j_tot_av, a_j_tot_av)
        )
        variables.update(self._get_standard_exchange_current_variables(j0))
        variables.update(self._get_standard_overpotential_variables(eta_r))

        variables.update(
            self._get_standard_volumetric_current_density_variables(variables)
        )

        # No SEI film resistance in this model
        eta_sei = pybamm.Scalar(0)
        variables.update(self._get_standard_sei_film_overpotential_variables(eta_sei))

        return variables

    def _get_diffusion_limited_current_density(self, variables):
        param = self.param
        if self.domain == "negative":
            if self.order == "leading":
                j_p = variables[
                    f"X-averaged positive electrode {self.reaction_name}"
                    "interfacial current density"
                ]
                j = -self.param.p.l * j_p / self.param.n.l
            elif self.order in ["composite", "full"]:
                tor_s = variables["Separator electrolyte transport efficiency"]
                c_ox_s = variables["Separator oxygen concentration"]
                N_ox_neg_sep_interface = (
                    -pybamm.boundary_value(tor_s, "left")
                    * param.curlyD_ox
                    * pybamm.BoundaryGradient(c_ox_s, "left")
                )
                N_ox_neg_sep_interface.domains = {"primary": "current collector"}

                j = -N_ox_neg_sep_interface / param.C_e / -param.s_ox_Ox / param.n.l

        return j
