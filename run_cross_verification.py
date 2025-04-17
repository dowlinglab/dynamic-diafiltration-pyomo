"""
Empirical B model verified with additional datasets
Xinhong Liu
University of Notre Dame
"""

from utility import *

# Parameter results from standard lag diafiltration experiment
# A.0
fit_stru_base = {'parameters': {'Lp': 11.113241068147595, 
            'beta_0': 1.0649294788103598, 
            'beta_1': 0.015171421224603474, 
            'sigma': 1.0}}

# Start cross verification
mode = 'DATA'
# A.1
data_stru = loadmat('data_library/data_stru-dataset270611.121.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# A.2
data_stru = loadmat('data_library/data_stru-dataset270711.121.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# B.1
data_stru = loadmat('data_library/data_stru-dataset270511.221.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# B.2
data_stru = loadmat('data_library/data_stru-dataset270511.321.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1, sigma_fixed=False)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# C.1
data_stru = loadmat('data_library/data_stru-dataset270511.421.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1, sigma_fixed=False)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# C.2
data_stru = loadmat('data_library/data_stru-dataset270511.921.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# D.1
data_stru = loadmat('data_library/data_stru-dataset270511.521.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# D.2
data_stru = loadmat('data_library/data_stru-dataset270511.621.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# E.1
data_stru = loadmat('data_library/data_stru-dataset270511.721.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1, sigma_fixed=False)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# E.2
data_stru = loadmat('data_library/data_stru-dataset270511.821.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)
