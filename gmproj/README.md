# Code for modeling battery degradation and expansion

Data files not included in this repo, please download from
 this [Google Drive Folder](https://drive.google.com/drive/folders/16uwOXhK_kvs6xNQBIiVQT5VzPDkkNnov?usp=sharing) and paste the files in the empty data folder provided. 
 
Ensure to paste the data in the corresponding subfolders.

- Run [run_model.ipynb](./run_model.ipynb) to simulate Aging for all cells at room temperature
- Run [fast_sim_exp](./fast_sim.ipynb) to simulate aging and expansion for C/5 cycling. Included differntial voltage and expansion analysis
- Run [predict_aging](./predict_aging.ipynb) to train and plot parameters using partial data window.
- Run [predict_aging_par](./predict_aging_par.ipynb) to train and plot parameters using partial data window with simulations running parallely to reduce simulation time. 
- Run [predict_aging_ECS_plots](./predict_aging_ECS_plots.ipynb) to generate plots for ECS 2022 including new analysis for predicing 50% DOD case 
