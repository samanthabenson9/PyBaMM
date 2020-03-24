#
# Base class for tab cooling submodels
#
import pybamm


class BaseModel(pybamm.BaseSubModel):
    """Base class for tab cooling submodels

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel


    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param):
        super().__init__(param)
