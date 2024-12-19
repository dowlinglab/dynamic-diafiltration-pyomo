"""
Diafiltration model variations and FIM calculations for DATA2 paper
Xinhong Liu
University of Notre Dame
"""

from utility import *

# Lag - w/ time correction <M1>
data_stru = loadmat('data_library/data_stru-dataset270511.123.mat')['data_stru']
mode = 'Lag'
theta = {'Lp': 11,
  'beta_c': 15,
  'beta_0': 1,
  'beta_1': 0.01,
  'sigma': 1.0,
  'S0': 0}
fit_stru_base1, sim_stru, sim_inter = solve_model(data_stru, mode, theta, sim_opt=False, B_form=1,LOUD=False)
file_name = os.path.join(mode+"M1-fit.json")
store_json(file_name,fit_stru_base1)

doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base1['parameters'], step=1e-8, formula='Backward', B_form=1)
file_name = os.path.join(mode+"M1-FIM.json")
store_json(file_name,doe_stru)

# Lag - w/o time correction <M2>
data_stru = loadmat('data_library/data_stru-dataset270511.122.mat')['data_stru']
mode = 'Lag'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, theta=fit_stru_base1['parameters'], sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M2-fit.json")
store_json(file_name,fit_stru)

# Lag truncated (DATA) - w/ time correction <M3>
data_stru = loadmat('data_library/data_stru-dataset270511.121.mat')['data_stru']
mode_truc = 'DATA'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode_truc, theta=fit_stru_base1['parameters'], sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M3-fit.json")
store_json(file_name,fit_stru)

doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base1['parameters'], step=1e-8, formula='Backward', B_form=1)
file_name = os.path.join(mode+"M3-FIM.json")
store_json(file_name,doe_stru)

# Lag truncated (DATA) - w/o time correction <M4>
data_stru = loadmat('data_library/data_stru-dataset270511.12.mat')['data_stru']
mode_truc = 'DATA'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode_truc, theta=fit_stru_base1['parameters'], sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M4-fit.json")
store_json(file_name,fit_stru)



# Overflow - w/ time correction <M1>
data_stru = loadmat('data_library/data_stru-dataset270511.423.mat')['data_stru']
mode = 'Overflow'
theta = {'Lp': 11,
  'beta_c': 15,
  'beta_0': 1,
  'beta_1': 0.01,
  'sigma': 1.0,
  'S0': -0.1756665334051245}
fit_stru_base2, sim_stru, sim_inter = solve_model(data_stru, mode, theta, sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M1-fit.json")
store_json(file_name,fit_stru_base2)

doe_stru = calc_FIM(data_stru, mode, theta=fit_stru_base2['parameters'], step=1e-8, formula='Backward', B_form=1)
file_name = os.path.join(mode+"M1-FIM.json")
store_json(file_name,doe_stru)
    
    
# Overflow - w/o time correction <M2>
data_stru = loadmat('data_library/data_stru-dataset270511.422.mat')['data_stru']
mode = 'Overflow'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, theta=fit_stru_base2['parameters'], sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M2-fit.json")
store_json(file_name,fit_stru)

# Overflow truncated (DATA) - w/ time correction <M3>
data_stru = loadmat('data_library/data_stru-dataset270511.421.mat')['data_stru']
mode_truc = 'DATA'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode_truc, theta=fit_stru_base2['parameters'], sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M3-fit.json")
store_json(file_name,fit_stru)

doe_stru = calc_FIM(data_stru, mode_truc, theta=fit_stru_base2['parameters'], step=1e-8, formula='Backward', B_form=1)
file_name = os.path.join(mode+"M3-FIM.json")
store_json(file_name,doe_stru)

# Overflow truncated (DATA) - w/o time correction <M4>
data_stru = loadmat('data_library/data_stru-dataset270511.42.mat')['data_stru']
mode_truc = 'DATA'
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode_truc, theta, sim_opt=False, B_form=1)
file_name = os.path.join(mode+"M4-fit.json")
store_json(file_name,fit_stru)
