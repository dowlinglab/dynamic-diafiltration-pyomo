# %% This code plots but does not process any data. This code uses data that has been been processed in MATLAB 
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib import patheffects
import pandas as pd

from diafiltration_plots import loadmat
from diafiltration_plots import plot_sim_comparison
from diafiltration_plots import plot_contour
from diafiltration_plots import plot_sim
from diafiltration_plots import plot_sim_show

# %% filtration parameter estimation(hybrid) using the data file 'data_stru-dataset501.1.mat'
data_stru_f_501_1 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.1.mat')['data_stru']
fit_stru_501_1 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.1 concpolar/fit_stru.mat')['fit_stru']
plot_sim_comparison(data_stru_f_501_1,fit_stru_501_1,plot_pred=True,cond=True,lg=False)

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

# # %% filtration parameter estimation(hybrid) using the data file 'data_stru-dataset501.11.mat'
# data_stru_f_501_11 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.11.mat')['data_stru']
# fit_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/fit_stru.mat')['fit_stru']
# plot_sim_comparison(data_stru_f_501_11,fit_stru,plot_pred=True,cond=True,lg=False)

# df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/contourdata-x_B-y_Lp.csv')
# plot_contour(df)
# plt.show()

# df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11 concpolar/contourdata-x_sigma-y_Lp.csv')
# plot_contour(df)
# plt.show()

# # filtration parameter estimation(hybrid)
# # data_stru_f_501_11 = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.11.mat')['data_stru']
# fit_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/fit_stru.mat')['fit_stru']
# plot_sim_comparison(data_stru_f_501_11,fit_stru,plot_pred=True,cond=False,lg=False)

# df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/contourdata-x_B-y_Lp.csv')
# plot_contour(df)
# plt.show()

# df = pd.read_csv(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/501.11/contourdata-x_sigma-y_Lp.csv')
# plot_contour(df)
# plt.show()

# # %% filtration sigma sensitivity
# sigma = [0.1,0.5,0.9]
# colorstring = 'rbg'

# sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.1.mat')['sim_stru']
# plot_sim(sim_stru,colorstring[0],'dashed')
# sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.5.mat')['sim_stru']
# plot_sim(sim_stru,colorstring[1],'solid')
# sim_stru = loadmat(r'/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/sigma sensitivity/sim_stru-dat501.1 C_Fin5.2843sig0.9.mat')['sim_stru']
# plot_sim(sim_stru,colorstring[2],'dotted')
# plot_sim_show(sigma,colorstring)