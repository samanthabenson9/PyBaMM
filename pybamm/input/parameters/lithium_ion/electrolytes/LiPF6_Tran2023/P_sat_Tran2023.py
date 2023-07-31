def P_sat_Tran2023(T):
    """
    Electrolyte saturation pressure as a function of temperature [K]

    References
    ----------
    .. [1] 

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte saturation pressure [kPa]
    """
    P_dmc = 10**(6.4338-1413.0/(T-44.25))*0.7
    P_ec = 10**(6.4897-1836.57/(T-102.23))*0.3
    P = P_ec + P_dmc 
    return P