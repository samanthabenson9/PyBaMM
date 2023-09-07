# Locating Model Parameters for updating

The model parameters are located in [input folder under pybamm ./pybamm/input/parameters/lithium_ion](../pybamm/input/parameters/lithium_ion/)

The folders for parameters are further divided into `cells`, `electrolytes`, `negative_electrodes`, `positive_electrodes` and `separators`

Under these folders, each parameter set has its own folder. The parameter set which we will be updating is called `UMBL_Siegel2022`

## Cell
Overall cell related parameters like nominal capacity, electrode thickness, current collector properties etc. are located in the cells folder [../pybamm/input/parameters/lithium_ion/cells/UMBL_Siegel2022/](../pybamm/input/parameters/lithium_ion/cells/UMBL_Siegel2022/)

Find [parameters.csv](../pybamm/input/parameters/lithium_ion/cells/UMBL_Siegel2022/parameters.csv) file in this folder and update the file with the new parameters.

## Electrolyte
Electrolyte related parameters are located in the electolytes folder [../pybamm/input/parameters/lithium_ion/electrolytes/UMBL_Siegel2022/](../pybamm/input/parameters/lithium_ion/electrolytes/LiPF6_Siegel2022/)

Find [parameters.csv](../pybamm/input/parameters/lithium_ion/electrolytes/LiPF6_Siegel2022/parameters.csv) file in this folder and update the file with the new parameters.

Additionally, the [conductivity](../pybamm/input/parameters/lithium_ion/electrolytes/LiPF6_Siegel2022/electrolyte_conductivity_Siegel.py) and [diffusivity](../pybamm/input/parameters/lithium_ion/electrolytes/LiPF6_Siegel2022/electrolyte_diffusivity_Siegel.py) functions are located in separate python files. Update them separately.

## Negative electrode
Negative electrode related parameters are located here [../pybamm/input/parameters/lithium_ion/negative_electrodes/graphite_UMBL_Siegel2022/](../pybamm/input/parameters/lithium_ion/negative_electrodes/graphite_UMBL_Siegel2022/)

Update these files in the folder
- `parameters.csv`
- `graphite_ocp_Siegel.py`
- `graphite_diffusivity_Siegel.py`
- `graphite_electrolyte_exchange_current_density_Siegel.py`
- `graphite_electrolyte_entropic_change_density_Siegel.py`

## Positive electrode
Positive electrode related parameters are located here [../pybamm/input/parameters/lithium_ion/positive_electrodes/NMC_UMBL_Siegel2022/](../pybamm/input/parameters/lithium_ion/positive_electrodes/NMC_UMBL_Siegel2022/)

Update these files in the folder
- `parameters.csv`
- `NMC_ocp_Siegel.py`
- `NMC_diffusivity_Siegel.py`
- `NMC_electrolyte_exchange_current_density_Siegel.py`
- `NMC_electrolyte_entropic_change_density_Siegel.py`

## Separator
Separator related parameters are located here [../pybamm/input/parameters/lithium_ion/separators/separator_Siegel2022/](../pybamm/input/parameters/lithium_ion/separators/separator_Siegel2022/)

Find [parameters.csv](../pybamm/input/parameters/lithium_ion/separators/separator_Siegel2022/parameters.csv) file in this folder and update the file with the new parameters.