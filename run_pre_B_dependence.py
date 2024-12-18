"""
Preliminary B dependence investigation
Xinhong Liu
University of Notre Dame
"""

from utility import *

mode = 'DATA'
Allsim_stru=dict()

data_stru = loadmat('data_library/data_stru-dataset270511.12.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['A.0'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270611.12.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['A.1'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270711.12.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['A.2'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.22.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['B.1'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.32.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['B.2'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.42.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['C.1'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.92.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['C.2'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.52.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['D.1'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.62.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['D.2'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.72.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['E.1'] = sim_stru

data_stru = loadmat('data_library/data_stru-dataset270511.82.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, sim_opt=False, B_form='pervial')
Allsim_stru['E.2'] = sim_stru


cmap1 = plt.get_cmap("tab10")
cmap2 = plt.get_cmap("tab20")
fig = plt.figure(figsize=(6,4))
fi = -1

for key, sim_st in Allsim_stru.items():
    if sim_st[0]['B'] != None:
        if fi == -1:
            co=cmap1(10)
        else:
            co=cmap2(fi)
        plt.plot([],[],color=co,linewidth=3,alpha=.8,label=key)
        for i in sim_st:
            if i==0:
                plt.plot(np.array(sim_st[i]['cIn'][80:]),np.array(sim_st[i]['Js'][80:])/np.array(sim_st[i]['Jw'][80:]),color=co,linewidth=2,
                      alpha=.8)
            else:
                plt.plot(np.array(sim_st[i]['cIn'][50:]),np.array(sim_st[i]['Js'][50:])/np.array(sim_st[i]['Jw'][50:]),color=co,linewidth=2,
                      alpha=.8)
        fi += 1

ylabelstr = 'Js/Jw [mM]'        
plt.xlabel('Interface Concentration[mM]',fontsize=16,fontweight='bold')
plt.ylabel(ylabelstr,fontsize=16,fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tick_params(direction="in")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend(fontsize=10,loc='best')

fname = 'Js_Jw_cin'
fig.savefig('figures/'+fname+'.png',dpi=300,bbox_inches='tight')

cmap1 = plt.get_cmap("tab10")
cmap2 = plt.get_cmap("tab20")
fig = plt.figure(figsize=(6,4))
fi = -1

for key, sim_st in Allsim_stru.items():
    if sim_st[0]['B'] != None:
        if fi == -1:
            co=cmap1(10)
        else:
            co=cmap2(fi)
        plt.plot([],[],color=co,linewidth=3,alpha=.8,label=key)
        for i in sim_st:    
            plt.plot(sim_st[i]['cIn'],sim_st[i]['B'],color=co,linewidth=2,
                      alpha=.8)
        fi += 1

ylabelstr = 'B [$\mathbf{\mu}$m $\mathbf{\cdot}$ s$\mathbf{^{-1}}$]'        
plt.xlabel('Interface Concentration[mM]',fontsize=16,fontweight='bold')
plt.ylabel(ylabelstr,fontsize=16,fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tick_params(direction="in")
plt.xlim(left=0)
plt.ylim(bottom=0)

fname = 'Bpervial'
fig.savefig('figures/'+fname+'.png',dpi=300,bbox_inches='tight')