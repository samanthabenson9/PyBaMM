#
# Class for graphite anode decomposition in Li-ion batteries
#
import pybamm
from scipy import constants


class PouchVenting(pybamm.BaseSubModel):
    """Pouch cell venting accounting for CO2 gas generation from SEI decomposition 
    and electrolyte saturation pressure.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    reactions : dict, optional
        Dictionary of reaction terms

    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param):
        super().__init__(param)

    def get_fundamental_variables(self):
        delta_sigma = pybamm.Variable("Cell expansion stress [kPa]", domain="current collector")

        variables = {
            "Cell expansion stress [kPa]": delta_sigma,
        }
        return variables

    def get_coupled_variables(self, variables):
        param = self.param
        m_an = param.therm.m_an
        x_sei_0 = param.therm.x_sei_0
        M_C6 = param.therm.M_an
        V_head_0 = param.therm.V_head_0
        A_surf = param.therm.A_surf
        L = param.therm.L_poron
        E = param.therm.E_poron
        alpha_cell = param.therm.alpha_cell
        P_crit = param.therm.P_crit
        sigma_0 = param.therm.sigma_0
        P_atm = param.therm.P_atm

        x_sei = variables["Fraction of Li in SEI"]
        delta_sigma = variables["Cell expansion stress [kPa]"]
        T = variables["X-averaged negative electrode temperature"]
        T_dimensional = param.Delta_T * T + param.T_ref
        delta_T = param.Delta_T * T

        delta_d = L*delta_sigma/E
        V_head = V_head_0 + A_surf*(delta_d-alpha_cell*delta_T)   
        P_sat = param.therm.P_sat(T) #dimensionless
        n_CO2 = m_an*(x_sei_0 - x_sei)/(2*M_C6)
        P_CO2 = (n_CO2*constants.R*T/V_head/1000)/P_crit
        P_total = P_CO2 + P_sat

        delta_sigma_gas = P_total-sigma_0-P_atm
        dleta_sigma_themal_expansion = E*alpha_cell*delta_T/L

        
        variables = {
            "Electrolyte gas saturation pressure [kPa]": P_sat*P_crit,
            "CO2 gas pressure [kPa]": P_CO2*P_crit,
            "Headspace volume [m3]": V_head,
            "Total gas pressure [kPa]": P_total*P_crit,
            
        }

        return variables

    def set_algebraic(self, variables):
        param = self.param
        P_crit = param.therm.P_crit
        L = param.therm.L_poron
        E = param.therm.E_poron
        alpha_cell = param.therm.alpha_cell
        P_atm = param.therm.P_atm
        sigma_0 = param.therm.sigma_0
    
        T = variables["X-averaged negative electrode temperature"]
        delta_T = param.Delta_T * T
        delta_sigma = variables["Cell expansion stress [kPa]"]
        P_total = variables["Total gas pressure [kPa]"]        
        self.algebraic = {delta_sigma: pybamm.maximum(P_total-sigma_0-P_atm, E*alpha_cell*delta_T/L)/P_crit} #pybamm.maximum(P_total-sigma_0-P_atm, E*alpha_cell*delta_T/L)

    def set_initial_conditions(self, variables):
        # param = self.param
        # sigma_0 = param.therm.sigma_0
        delta_sigma = variables["Cell expansion stress [kPa]"]
        self.initial_conditions = {
            delta_sigma: pybamm.Scalar(0), 
            }
        
    def set_events(self, variables):
        """
        A method to set events related to the state of submodel variable. Note: this
        method modifies the state of self.events. Unless overwritten by a submodel, the
        default behaviour of 'pass' is used a implemented in
        :class:`pybamm.BaseSubModel`.

        Parameters
        ----------
        variables: dict
            The variables in the whole model.
        """

        delta_sigma = variables["Cell expansion stress [kPa]"]
        self.events.append(
            pybamm.Event(
                "Switch expansion mode",
                (P_total-sigma_0-P_atm) - (E*alpha_cell*delta_T/L),
                pybamm.EventType.DISCONTINUITY,
            )
        )
