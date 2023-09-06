# Code for modeling battery degradation and expansion

## Installation instructions
1. Install [git](https://git-scm.com/downloads) 
2. Install [python 3.9](https://www.python.org/downloads/release/python-3913/)
3. Clone the git repo:
```
git clone https://github.com/js1tr3/PyBaMM
```
4. Follow install from source instruction from PyBaMM documentation to install our custom fork on PyBaMM
https://docs.pybamm.org/en/latest/source/user_guide/installation/GNU-linux.html
5. Switch to the `deg-model-pub` branch of the PyBaMM fork
```
git checkout deg-model-pub
```

6. Change the working directory to `degradation_model` subfolder
```
cd degradation_model
```
7. Data files not included in this repo, please download from this [Google Drive Folder](https://drive.google.com/drive/folders/16uwOXhK_kvs6xNQBIiVQT5VzPDkkNnov?usp=sharing) and paste the files in the empty data folder provided. Ensure to paste the data in the corresponding subfolders.
8.  Run [run_model.ipynb](./run_model.ipynb) notebook to simulate Aging for all cells at room temperature
  - Includes Resistance simulations
  - Includes Expansion simulations
  - Generates figures in the results section of the paper. 
9.  Run [other_figures.ipynb](./run_model.ipynb) to generate figures from other sections in the paper