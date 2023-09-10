import pybamm


def NMC_ocp_AndrewMPM(x):
    """
    Nickel Managanese Cobalt Oxide (NMC) Open Circuit Potential (OCP) as a
    function of the stochiometry. The fit is taken from Peyman MPM.

    References
    ----------
    Peyman MPM manuscript (to be submitted)

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
       Stochiometry of material (li-fraction)

    """
# Peyman
    
#     u_eq = (
#         4.3452
#         - 1.6518 * sto
#         + 1.6225 * (sto ** 2)
#         - 2.0843 * (sto ** 3)
#         + 3.5146 * (sto ** 4)
#         - 2.2166 * (sto ** 5)
#         - 0.5623 * pybamm.exp(109.451 * sto - 100.006)
#     )
    
#     Andrew_not shifted
#     u_eq = (
#         2.992
#         - 2.098 * sto
#         - 0.6943 * (sto ** 2)
#         + 4.341 * (sto ** 3)
#         - 3.883 * (sto ** 4)
#         + 0.611 * (sto ** 5)
#         + 0.8258 * pybamm.exp(0.4484 * sto + 0.4757)
#     )
    
#   GM July 2022
    p1 =  -199.3069
    p2 =  762.6359
    p3 =  -1.1957e+03
    p4 =   1.0031e+03
    p5 =  -500.7007
    p6 =  158.2556
    p7 = -32.0707
    p8 =   4.6323
    p9 =  -1.6569
    p10 =  4.2932

    u_eq = p1*x**9 + p2*x**8 + p3*x**7 + p4*x**6 + p5*x**5 + p6*x**4 + p7*x**3 + p8*x**2 + p9*x + p10;    

    return u_eq


# if __name__ == "__main__":  # pragma: no cover
#     x = pybamm.linspace(0, 1)
#     pybamm.plot(x, NMC_ocp_PeymanMPM(x))
