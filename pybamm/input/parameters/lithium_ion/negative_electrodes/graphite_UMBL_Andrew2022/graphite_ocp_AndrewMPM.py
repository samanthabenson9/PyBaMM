import pybamm


def graphite_ocp_AndrewMPM(sto):
    """
    Graphite Open Circuit Potential (OCP) as a function of the
    stochiometry. The fit is taken from Peyman MPM [1].

    References
    ----------
    .. [1] Peyman Mohtat et al, MPM (to be submitted)
    """
# Peyman
    u_eq = (
#         0.063
#         + 0.8 * pybamm.exp(-75 * (sto + 0.007))
#         - 0.0120 * pybamm.tanh((sto - 0.127) / 0.016)
#         - 0.0118 * pybamm.tanh((sto - 0.155) / 0.016)
#         - 0.0035 * pybamm.tanh((sto - 0.220) / 0.020)
#         - 0.0095 * pybamm.tanh((sto - 0.190) / 0.013)
#         - 0.0145 * pybamm.tanh((sto - 0.490) / 0.020)
#         - 0.0800 * pybamm.tanh((sto - 1.030) / 0.055)
#  Andrew       
        0.09753
        + 0.6722 * pybamm.exp(-41.1 * (sto + 0.0029))
        - 0.0147 * pybamm.tanh((sto - 0.1503) / 0.0161)
        - 0.0121 * pybamm.tanh((sto - 0.1789) / 0.0179)
        - 0.0024 * pybamm.tanh((sto - 0.2503) / 0.0069)
        - 0.0083 * pybamm.tanh((sto - 0.2118) / 0.0133)
        - 0.0125 * pybamm.tanh((sto - 0.5158) / 0.0180)
        - 0.0538 * pybamm.tanh((sto - 0.9961) / 0.0482)
        
        
    )

    return u_eq


# if __name__ == "__main__":  # pragma: no cover
#     x = pybamm.linspace(1e-10, 1 - 1e-10, 1000)
#     # pybamm.plot(x, graphite_ocp_PeymanMPM(x))
#     pybamm.plot(x, -1e-8 * pybamm.log(x / (1 - x)))
