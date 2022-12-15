import warnings
import importlib_metadata
import textwrap
from collections.abc import Mapping


class ParameterSets(Mapping):
    """
    Dict-like interface for accessing registered pybamm parameter sets.
    Access via :py:data:`pybamm.parameter_sets`

    Examples
    --------
    Listing available parameter sets:

    .. doctest::

        >>> import pybamm
        >>> list(pybamm.parameter_sets)
        ['Ai2020', 'Chen2020', ...]

    Get the docstring for a parameter set:

    .. doctest::

        >>> import pybamm
        >>> print(pybamm.parameter_sets.get_docstring("Ai2020"))
        <BLANKLINE>
        Parameters for the Enertech cell (Ai2020), from the papers:
        ...

<<<<<<< HEAD
Mohtat2020 = {
    "chemistry": "lithium_ion",
    "cell": "UMBL_Mohtat2020",
    "negative electrode": "graphite_UMBL_Mohtat2020",
    "separator": "separator_Mohtat2020",
    "positive electrode": "NMC_UMBL_Mohtat2020",
    "electrolyte": "LiPF6_Mohtat2020",
    "experiment": "1C_charge_from_empty_Mohtat2020",
    "sei": "example",
    "lithium plating": "mohtat2020_Li_plating",
    "citation": "Mohtat2020",
}

Siegel2022 = {
    "chemistry": "lithium_ion",
    "cell": "UMBL_Siegel2022",
    "negative electrode": "graphite_UMBL_Siegel2022",
    "separator": "separator_Siegel2022",
    "positive electrode": "NMC_UMBL_Siegel2022",
    "electrolyte": "LiPF6_Siegel2022",
    "experiment": "1C_charge_from_empty_Mohtat2020",
    "sei": "example",
    "lithium plating": "yang2017_Li_plating",
    # "citation": "Siegel2022",
}

Ramadass2004 = {
    "chemistry": "lithium_ion",
    "cell": "sony_Ramadass2004",
    "negative electrode": "graphite_Ramadass2004",
    "separator": "separator_Ecker2015",  # no values found, relevance?
    "positive electrode": "lico2_Ramadass2004",
    "electrolyte": "lipf6_Ramadass2004",
    "experiment": "1C_discharge_from_full_Ramadass2004",
    "sei": "ramadass2004",
    "citation": "Ramadass2004",
}
=======
    See also: :ref:`adding-parameter-sets`

    """
>>>>>>> upstream/develop

    def __init__(self):
        # Dict of entry points for parameter sets, lazily load entry points as
        self.__all_parameter_sets = dict()
        for entry_point in importlib_metadata.entry_points(
            group="pybamm_parameter_sets"
        ):
            self.__all_parameter_sets[entry_point.name] = entry_point

    def __new__(cls):
        """Ensure only one instance of ParameterSets exists"""
        if not hasattr(cls, "instance"):
            cls.instance = super(ParameterSets, cls).__new__(cls)
        return cls.instance

    def __getitem__(self, key) -> dict:
        return self.__load_entry_point__(key)()

    def __load_entry_point__(self, key) -> callable:
        """Check that ``key`` is a registered ``pybamm_parameter_sets``,
        and return the entry point for the parameter set, loading it needed.
        """
        if key not in self.__all_parameter_sets:
            raise KeyError(f"Unknown parameter set: {key}")
        ps = self.__all_parameter_sets[key]
        if isinstance(ps, importlib_metadata.EntryPoint):
            ps = self.__all_parameter_sets[key] = ps.load()
        return ps

    def __iter__(self):
        return self.__all_parameter_sets.__iter__()

    def __len__(self) -> int:
        return len(self.__all_parameter_sets)

    def get_docstring(self, key):
        """Return the docstring for the ``key`` parameter set"""
        return textwrap.dedent(self.__load_entry_point__(key).__doc__)

    def __getattribute__(self, name):
        try:
            return super().__getattribute__(name)
        except AttributeError as error:
            # For backwards compatibility, parameter sets that used to be defined in
            # this file now return the name as a string, which will load the same
            # parameter set as before when passed to `ParameterValues`
            if name in self:
                msg = (
                    "Parameter sets should be called directly by their name ({0}), "
                    "instead of via pybamm.parameter_sets (pybamm.parameter_sets.{0})."
                ).format(name)
                warnings.warn(msg, DeprecationWarning)
                return name
            raise error


#: Singleton Instance of :class:ParameterSets """
parameter_sets = ParameterSets()
