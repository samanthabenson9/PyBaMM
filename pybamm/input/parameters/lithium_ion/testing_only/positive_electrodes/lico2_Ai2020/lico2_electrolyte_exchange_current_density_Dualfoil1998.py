from pybamm import exp, constants, Parameter


def lico2_electrolyte_exchange_current_density_Dualfoil1998(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between lico2 and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [2] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
<<<<<<< HEAD:pybamm/input/parameters/lithium_ion/negative_electrodes/graphite_UMBL_Mohtat2020/graphite_electrolyte_exchange_current_density_PeymanMPM.py
    m_ref =  Parameter("Negative electrode reference exchange-current density [A.m-2(m3.mol)1.5]")
    # m_ref = 4*1.061 * 10 ** (-6)  # unit has been converted
    # units are (A/m2)(mol/m3)**1.5 - includes ref concentrations
    E_r = 37480
=======
    m_ref = 1 * 10 ** (-11) * constants.F  # need to match the unit from m/s
    # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    E_r = 5000
>>>>>>> upstream/develop:pybamm/input/parameters/lithium_ion/testing_only/positive_electrodes/lico2_Ai2020/lico2_electrolyte_exchange_current_density_Dualfoil1998.py
    arrhenius = exp(E_r / constants.R * (1 / 298.15 - 1 / T))

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )
