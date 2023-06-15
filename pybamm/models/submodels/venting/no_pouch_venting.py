#
# Class for when there is no anode decomposition
#
import pybamm
from scipy import constants


class NoPouchVenting(pybamm.BaseSubModel):
    """Base class for no pouch cell venting.

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

        variables = {
            "Electrolyte gas saturation pressure [kPa]": pybamm.Scalar(0),
            "CO2 gas pressure [kPa]": pybamm.Scalar(0),
            "Headspace volume [m3]": pybamm.Scalar(0),
            "Expansion stress [kPa]": pybamm.Scalar(0),
        }

        return variables
