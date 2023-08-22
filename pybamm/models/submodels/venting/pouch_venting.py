#
# Class for graphite anode decomposition in Li-ion batteries
#
import pybamm
from scipy import constants


class PouchVenting(pybamm.BaseSubModel):
    """Pouch cell venting accounting for CO2 gas generation from SEI decomposition 
    and electrolyte saturation pressure. Applies electrolyte dryout based on Kupper et al 2018. 

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
        # a = pybamm.Variable("Effective loss of active material") # Electrode activity
        # s = pybamm.Variable("Electrolyte saturation") # Liquid saturation 
        # eps_p_s = pybamm.Variable("Positive electrode active material volume fraction")
        # eps_n_s = pybamm.Variable("Negative electrode active material volume fraction")
        variables = { 
            "Cell expansion stress [kPa]": delta_sigma,
            # "Electrolyte saturation": s,
            # "Effective loss of active material": a,
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
        M_e = param.therm.M_e
        rho_e = param.therm.rho_e
        m_e_0 = param.therm.m_e_0
        epsilon_p_s = param.p.prim.epsilon_s
        epsilon_n_s = param.n.prim.epsilon_s

        x_sei = variables["Fraction of Li in SEI"]
        delta_sigma = variables["Cell expansion stress [kPa]"]
        # a = variables["Effective loss of active material"]
        T = variables["X-averaged negative electrode temperature"] 
        T_dimensional = param.Delta_T * T + param.T_ref
        # T_amb_dim = param.T_amb_dim(pybamm.t * param.timescale)
        delta_T = param.Delta_T * T

        delta_d = L/E*delta_sigma
        V_head = V_head_0 + A_surf*(delta_d-alpha_cell*delta_T)   
        P_sat = param.therm.P_sat_dimensional(T_dimensional) 
        P_sat_ref = param.therm.P_sat_dimensional(param.T_ref) 

        n_CO2 = m_an*(x_sei_0 - x_sei)/(2*M_C6)
        P_CO2 = (n_CO2*constants.R*T_dimensional/V_head/1000) #kPa
        P_total = P_CO2 + P_sat 

        # delta_sigma_gas = P_total-sigma_0-P_atm
        # delta_sigma_themal_expansion = E*alpha_cell*delta_T/L

        n_e_vap = P_sat*1000*V_head/(constants.R*T_dimensional) - P_sat_ref*1000*V_head_0/(constants.R*param.T_ref)
        m_e_loss = pybamm.maximum(n_e_vap*M_e,0)  # [kg]
        V_e_loss = m_e_loss/rho_e  # [m3]

        s_loss = m_e_loss/m_e_0
        s =  (1-s_loss*400)
        a = ((0.5 * s + 0.5) * (0.5 * pybamm.tanh(15 * (s - 0.4)) + 0.5))

        variables.update({
            "Electrolyte gas saturation pressure": P_sat/P_crit,
            "CO2 gas pressure": P_CO2/P_crit,
            "Total gas pressure": P_total/P_crit,
            "Headspace volume": V_head/V_head_0,
            "Loss of electrolyte": m_e_loss/m_e_0,
            # "X-averaged positive electrode active material volume fraction": epsilon_p_s*a,
            # "X-averaged negative electrode active material volume fraction": epsilon_n_s*a,

            "Electrolyte gas saturation pressure [kPa]": P_sat,
            "CO2 gas pressure [kPa]": P_CO2,
            "Total gas pressure [kPa]": P_total,
            "Headspace volume [m3]": V_head,
            "Mass of electrolyte lost [kg]": m_e_loss,
            "Volume of electrolyte lost [m3]": V_e_loss,
        })

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
        # a = variables["Effective loss of active material"]
        # s = variables["Electrolyte saturation"]
        # s_loss = variables["Loss of electrolyte"]   
        self.algebraic = {
            delta_sigma: delta_sigma - pybamm.maximum(P_total-sigma_0-P_atm, E*alpha_cell*delta_T/L),
            # s: s-(1-s_loss*200), 
            # a: a- ((0.5 * s + 0.5) * (0.5 * pybamm.tanh(15 * (s - 0.4)) + 0.5))
            } 
        
    def set_initial_conditions(self, variables):
        # param = self.param
        # sigma_0 = param.therm.sigma_0
        delta_sigma = variables["Cell expansion stress [kPa]"]
        # a = variables["Effective loss of active material"]
        # s = variables["Electrolyte saturation"] 
        
        self.initial_conditions = {
            delta_sigma: pybamm.Scalar(0), 
            # a: pybamm.Scalar(1),
            # s: pybamm.Scalar(1),
            }
        
    # def set_events(self, variables):
    #     """
    #     A method to set events related to the state of submodel variable. Note: this
    #     method modifies the state of self.events. Unless overwritten by a submodel, the
    #     default behaviour of 'pass' is used a implemented in
    #     :class:`pybamm.BaseSubModel`.

    #     Parameters
    #     ----------
    #     variables: dict
    #         The variables in the whole model.
    #     """
    #     param = self.param
    #     sigma_0 = param.therm.sigma_0
    #     P_atm = param.therm.P_atm
    #     E = param.therm.E_poron
    #     L = param.therm.L_poron
    #     alpha_cell = param.therm.alpha_cell
    #     T = variables["X-averaged negative electrode temperature"]
    #     delta_T = param.Delta_T * T
    #     # delta_sigma = variables["Cell expansion stress [kPa]"]
    #     P_total = variables["Total gas pressure [kPa]"]
    #     self.events.append(
    #         pybamm.Event(
    #             "Switch expansion mode",
    #             (P_total-sigma_0-P_atm) - (E*alpha_cell*delta_T/L),
    #             pybamm.EventType.DISCONTINUITY,
    #         )
    #     )
