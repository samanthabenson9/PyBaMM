#
# Class for varying active material volume fraction, driven by stress
#
import pybamm

from .base_active_material import BaseModel


class ElectrolyteDryout(BaseModel):
    """Submodel for varying active material volume fraction from [1]_ and [2]_.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    options : dict
        Additional options to pass to the model
    x_average : bool
        Whether to use x-averaged variables (SPM, SPMe, etc) or full variables (DFN)

    **Extends:** :class:`pybamm.active_material.BaseModel`

    References
    ----------
    .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2019). Electrochemical
           Thermal-Mechanical Modelling of Stress Inhomogeneity in Lithium-Ion Pouch
           Cells. Journal of The Electrochemical Society, 167(1), 013512.
    .. [2] Reniers, J. M., Mulder, G., & Howey, D. A. (2019). Review and performance
           comparison of mechanical-chemical degradation models for lithium-ion
           batteries. Journal of The Electrochemical Society, 166(14), A3189.
    """

    def __init__(self, param, domain, options, x_average):
        super().__init__(param, domain, options=options)
        pybamm.citations.register("Reniers2019")
        self.x_average = x_average

    def get_fundamental_variables(self):
        domain = self.domain.lower() + " electrode"
        if self.x_average is True:
            eps_solid_xav = pybamm.Variable(
                "X-averaged " + domain + " active material volume fraction",
                domain="current collector",
            )
            eps_solid = pybamm.PrimaryBroadcast(eps_solid_xav, domain)
        else:
            eps_solid = pybamm.Variable(
                self.domain + " electrode active material volume fraction",
                domain=domain,
                auxiliary_domains={"secondary": "current collector"},
            )
        variables = self._get_standard_active_material_variables(eps_solid)

        # a = pybamm.Variable("Effective loss of active material") # Electrode activity
        # s = pybamm.Variable("Electrolyte saturation") # Liquid saturation 
        # variables.update({            
            # "Electrolyte saturation": s,
            # "Effective loss of active material": a,
            # })
        return variables

    def get_coupled_variables(self, variables):
        return variables


    def set_algebraic(self, variables):
        eps_solid_init = self.domain_param.prim.epsilon_s
        s_loss = variables["Loss of electrolyte"]
        s = (1-s_loss*100) #(1-s_loss*200)
        a0 = ((0.5 * s + 0.5) * (0.5 * pybamm.tanh(15 * (s - 0.4)) + 0.5))
        Domain = self.domain + " electrode"
        if self.x_average is True:
            eps_solid = variables[
                "X-averaged " + Domain.lower() + " active material volume fraction"
            ]
            eps_solid_init = pybamm.x_average(eps_solid_init)
            a = a0
        else:
            eps_solid = variables[Domain + " active material volume fraction"]
            a = pybamm.PrimaryBroadcast(a0, Domain.lower())

        self.algebraic = { 
        eps_solid: eps_solid - a*eps_solid_init
        } 

    def set_initial_conditions(self, variables):

        eps_solid_init = self.domain_param.prim.epsilon_s
        # a = variables["Effective loss of active material"]
        # s = variables["Electrolyte saturation"] 
        if self.x_average is True:
            eps_solid_xav = variables[
                "X-averaged "
                + self.domain.lower()
                + " electrode active material volume fraction"
            ]
            self.initial_conditions = {
                eps_solid_xav: pybamm.x_average(eps_solid_init),
                # a: pybamm.Scalar(1),
                # s: pybamm.Scalar(1),
            }
        else:
            eps_solid = variables[
                self.domain + " electrode active material volume fraction"
            ]
            self.initial_conditions = {
                eps_solid: eps_solid_init,
                # a: pybamm.Scalar(1),
                # s: pybamm.Scalar(1),
            }
