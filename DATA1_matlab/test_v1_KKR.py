# %% This code plots but does not process any data. This code uses data that has been been processed in MATLAB 

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib import patheffects
import pandas as pd

from diafiltration_plots import loadmat # function that loads experimental and simulated data
from diafiltration_plots import plot_sim_comparison # function to generate Fig. 2
from diafiltration_plots import plot_contour
from diafiltration_plots import plot_sim
from diafiltration_plots import plot_sim_show
from diafiltration_plots import plot_conc_range # function to generate Fig. 3

# %% To run the files associated with this code, either run the code as written
#    if you are certain that the directory/folder that the code is stored and 
#    run in is the same as where the 'data' folder is stored.
#    If you are uncertain of this or if teh code resides in a different location
#    provide the full path to the 'data' folder location.
#    For example: r'/Users/user_name/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/data_stru-dataset501.1.mat'


# %% Figure 2: Diafiltration (dataset used --> 'data_stru-dataset511.12.mat') 
#    and filtration (dataset used --> 'data_stru-dataset501.1.mat') mass and concentration profile

#    Filtration
# data_stru_f_501_1 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/data_stru-dataset501.1.mat')['data_stru']
data_stru_f_501_1 = loadmat(r'./data/data_stru-dataset501.1.mat')['data_stru']
# fit_stru_501_1 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/501.1 concpolar/fit_stru.mat')['fit_stru']
fit_stru_f_501_1 = loadmat(r'./data/501.1 concpolar/fit_stru.mat')['fit_stru']

plot_sim_comparison(data_stru_f_501_1, fit_stru_f_501_1, plot_pred = True, cond = True, lg = False)

#    Diafiltration
# data_stru_d_511_12 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset511.12.mat')['data_stru']
data_stru_d_511_12 = loadmat(r'./data/data_stru-dataset511.12.mat')['data_stru']
# fit_stru_d_511_12 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/511.12 concpolar/fit_stru.mat')['fit_stru']
fit_stru_d_511_12 = loadmat(r'./data/511.12 concpolar/fit_stru.mat')['fit_stru']

plot_sim_comparison(data_stru_d_511_12, fit_stru_d_511_12, plot_pred = True, cond = True, lg = False)

# %% Figure 3: Diafiltration (dataset used --> 'diafiltration.csv') 
#    and filtration (dataset used --> 'filtration.csv') phase space 
#    retentate vs. permeate concentration profile

# df_fil_ret_vs_per_conc = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/experiment space/filtration.csv',header=2)
df_fil_ret_vs_per_conc = pd.read_csv(r'./data/experiment space/filtration.csv',header=2)
# df_diafil_ret_vs_per_conc = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/experiment space/diafiltration.csv',header=2)
df_diafil_ret_vs_per_conc = pd.read_csv(r'./data/experiment space/diafiltration.csv',header=2)

plot_conc_range(df_fil_ret_vs_per_conc, df_diafil_ret_vs_per_conc) 

# %% Figure 4: Sigma sensitivity analysis with membrane M1 (501.1)

# filtration sigma sensitivity
# plt.figure()

sigma_f = [0.1,0.5,0.9]
colorstring_f = 'rbg'

# sim_stru_f_sigma_sens_501_1_sig_0_1 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.1.mat')['sim_stru']
sim_stru_f_sigma_sens_501_1_sig_0_1 = loadmat(r'./data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.1.mat')['sim_stru']
plot_sim(sim_stru_f_sigma_sens_501_1_sig_0_1, colorstring_f[0],'dashed')

# plt.figure()

# sim_stru_f_sigma_sens_501_1_sig_0_5 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.5.mat')['sim_stru']
sim_stru_f_sigma_sens_501_1_sig_0_5 = loadmat(r'./data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.5.mat')['sim_stru']
plot_sim(sim_stru_f_sigma_sens_501_1_sig_0_5, colorstring_f[1],'solid')

# plt.figure()

# sim_stru_f_sigma_sens_501_1_sig_0_9 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.9.mat')['sim_stru']
sim_stru_f_sigma_sens_501_1_sig_0_9 = loadmat(r'./data/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.9.mat')['sim_stru']
plot_sim(sim_stru_f_sigma_sens_501_1_sig_0_9, colorstring_f[2],'dotted')

plot_sim_show(sigma_f,colorstring_f)

# plt.figure()

# diafiltration sigma sensitivity
sigma_d = [0.1,0.5,0.9]
colorstring_d = 'rbg'

sim_stru_d_sigma_sens_511_12_sig_0_1 = loadmat(r'./data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.1.mat')['sim_stru']
plot_sim(sim_stru_d_sigma_sens_511_12_sig_0_1, colorstring_d[0],'dashed')

# plt.figure()

sim_stru_d_sigma_sens_511_12_sig_0_5  = loadmat(r'./data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.5.mat')['sim_stru']
plot_sim(sim_stru_d_sigma_sens_511_12_sig_0_5, colorstring_d[1],'solid')

# plt.figure()

sim_stru_d_sigma_sens_511_12_sig_0_9 = loadmat(r'./data/sigma sensitivity/sim_stru-dat511.12 C_Fin15.2052sig0.9.mat')['sim_stru']
plot_sim(sim_stru_d_sigma_sens_511_12_sig_0_9, colorstring_d[2],'dotted')

plot_sim_show(sigma_d, colorstring_d)

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.1 concpolar/contourdata-x_B-y_Lp.csv')
plot_contour(df)
plt.show()

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.1 concpolar/contourdata-x_sigma-y_Lp.csv')
plot_contour(df)
plt.show()

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.1/contourdata-x_B-y_Lp.csv')
plot_contour(df)
plt.show()

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.1/contourdata-x_sigma-y_Lp.csv')
plot_contour(df)
plt.show()

# %% filtration parameter estimation(hybrid) using the data file 'data_stru-dataset501.11.mat'
data_stru_f_501_11 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.11.mat')['data_stru']
fit_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/fit_stru.mat')['fit_stru']
plot_sim_comparison(data_stru_f_501_11,fit_stru,plot_pred=True,cond=True,lg=False)

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/contourdata-x_B-y_Lp.csv')
plot_contour(df)
plt.show()

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/contourdata-x_sigma-y_Lp.csv')
plot_contour(df)
plt.show()

# filtration parameter estimation(hybrid)
# data_stru_f_501_11 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.11.mat')['data_stru']
fit_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/fit_stru.mat')['fit_stru']
plot_sim_comparison(data_stru_f_501_11,fit_stru,plot_pred=True,cond=False,lg=False)

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/contourdata-x_B-y_Lp.csv')
plot_contour(df)
plt.show()

df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/contourdata-x_sigma-y_Lp.csv')
plot_contour(df)
plt.show()

# %% filtration sigma sensitivity
sigma = [0.1,0.5,0.9]
colorstring = 'rbg'

sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.1.mat')['sim_stru']
plot_sim(sim_stru,colorstring[0],'dashed')
sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.5.mat')['sim_stru']
plot_sim(sim_stru,colorstring[1],'solid')
sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.9.mat')['sim_stru']
plot_sim(sim_stru,colorstring[2],'dotted')
plot_sim_show(sigma,colorstring)