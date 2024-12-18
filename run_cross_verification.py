"""
Empirical B model verified with additional datasets
Xinhong Liu
University of Notre Dame
"""

from utility import *

def solve_model_B_fix(data_stru, mode, theta=None, sim_opt=False, B_form=1, LOUD=False):
    """
    Solve pyomo model
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        sim_opt: boolean, if run simulation with fixed parameter values
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration 
        LOUD: boolean, if print out parameter results and store model predictions
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
        sim_inter: dict, model predictions for experimental measurements 
    """
    
    # interpolation function
    def interpolation(m,var,n_vial,t):
        """
        Interpolation function for pyomo model
        """
        tp = list(m.tau)
        for j in range(0,len(tp)-1):
            if tp[j+1]>=t and tp[j]<=t:
                var_inter = (var[n_vial,tp[j+1]]-var[n_vial,tp[j]]) * (t - tp[j])/(tp[j+1]-tp[j]) + var[n_vial,tp[j]]
        return var_inter
    
    # objective function
    def obj_rule(m):
        """
        pyomo Objective function (residual calculated using interpolation from model)
        """
        obj_m = 0
        obj_cp = 0
        obj_cf0 = 0
        obj_cf = 0
        ob_m = 0
        ob_cp = 0
        ob_cf0 = 0
        ob_cf = 0
        collect_vial = data_stru['data_config']['n'] - data_stru['data_config']['n_extra']
        Count_m = 0
        Count_cp = 0
        Count_cf = 0
        count_cf0 = 0
        t_delay = data_stru['data_raw'][0]['time'][0]
        TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(data_stru['data_config']['n'])]
        TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
        TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(data_stru['data_config']['n'])]
        TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    
        for n_vial in m.n_vial:
            t_meas = data_stru['data_raw'][n_vial-1]['time'] - t_delay
            mv_meas = data_stru['data_raw'][n_vial-1]['mass']
            cp_meas = data_stru['data_raw'][n_vial-1]['cV_avg']
            cf_meas = data_stru['data_raw'][n_vial-1]['cF_exp']
    
            t_meas_scaled = [(t-TI_dict[n_vial])/(TF_dict[n_vial]-TI_dict[n_vial]) for t in t_meas]
            
            mv_pred=[]
            obj_mi = 0
            ob_mi = 0
            count_m = 0
            for i in range(0,len(t_meas_scaled)):
                mv_inter = interpolation(m,m.mV,n_vial,t_meas_scaled[i])
                mv_pred.append(mv_inter)
                
                if not np.isnan(mv_meas[i]):
                    res_m = mv_pred[i]-mv_meas[i]             
                    count_m += 1
                    obj_mi += (res_m/0.01)**2 # 0.01g error
                    ob_mi += res_m**2
            if n_vial >= data_stru['data_config']['n_v0']: 
                Count_m += count_m
                obj_m += obj_mi#/count_m/collect_vial
                ob_m += ob_mi
    
            # permeate residual squared / uncertainty / # of measurements
            if type(data_stru['data_raw'][n_vial-1]['cV_avg']) != list:
                if n_vial > data_stru['data_config']['n_extra']:
                    res_cp = m.cV[n_vial,m.tau.last()]-data_stru['data_raw'][n_vial-1]['cV_avg']
                    obj_cp += (res_cp/(0.03*data_stru['data_raw'][n_vial-1]['cV_avg']))**2 # 3% error
                    ob_cp += (res_cp/(data_stru['data_raw'][n_vial-1]['cV_avg']))**2
                    Count_cp = collect_vial
            else:
                cp_pred=[]
                obj_cpi = 0
                ob_cpi = 0
                count_cp = 0
                if not all(np.isnan(cp_meas)):
                    for i in range(0,len(t_meas_scaled)):
                        cp_inter = interpolation(m,m.cV,n_vial,t_meas_scaled[i])
                        cp_pred.append(cp_inter)
                        if not np.isnan(cp_meas[i]):   
                            res_cp = cp_pred[i]-cp_meas[i]
                            count_cp += 1
                            obj_cpi += (res_cp/(0.03*cp_meas[i]))**2 # 3% error
                            ob_cpi += (res_cp/(cp_meas[i]))**2
                    Count_cp = collect_vial
                    obj_cp += obj_cpi#/count_cp/collect_vial
                    ob_cp += ob_cpi
         
            # retentate residual squared / uncertainty / # of measurements
            cf_pred=[]
            obj_cfi = 0
            ob_cfi = 0
            count_cf = 0
            if not all(np.isnan(cf_meas)):
                for i in range(0,len(t_meas_scaled)):
                    cf_inter = interpolation(m,m.cF,n_vial,t_meas_scaled[i])
                    cf_pred.append(cf_inter)
                    if not np.isnan(cf_meas[i]):   
                        res_cf = cf_pred[i]-cf_meas[i]
                        count_cf += 1
                        obj_cfi += (res_cf/(0.003*cf_meas[i]))**2 # 0.3% error
                        ob_cfi += (res_cf/(cf_meas[i]))**2
                if n_vial >= data_stru['data_config']['n_v0']:
                    Count_cf += count_cf
                    obj_cf += obj_cfi#/count_cf/collect_vial
                    ob_cf += ob_cfi
                else: 
                    count_cf0 += count_cf
                    obj_cf0 += obj_cfi#/count_cf/collect_vial
                    ob_cf0 += ob_cfi
        
        m.count_m = Count_m
        m.count_cv = Count_cp
        m.count_cr = Count_cf
        m.count_cr0 = count_cf0
        m.count = Count_m+Count_cp+Count_cf+count_cf0
        
        m.obj_m = obj_m/Count_m
        m.obj_cv = obj_cp/Count_cp
        m.obj_cr = (obj_cf0+obj_cf)/(count_cf0+Count_cf)
        m.obj_cr_tru = obj_cf/Count_cf
        m.obj_tru = 1e4*(obj_m/Count_m+obj_cp/Count_cp+obj_cf/Count_cf)
        m.llh1 = m.count_m*log(m.obj_m) + m.count_cv*log(m.obj_cv) + (m.count_cr0+m.count_cr)*log(m.obj_cr)#m.count * log((obj_m+obj_cp+obj_cf0+obj_cf)/m.count)
        m.llh2 = Count_m*log(ob_m/Count_m) + Count_cp*log(ob_cp/Count_cp) + (count_cf0+Count_cf)*log((ob_cf0+ob_cf)/(count_cf0+Count_cf))
        
        return 1e4*(obj_m/Count_m + obj_cp/Count_cp + (obj_cf0+obj_cf)/(count_cf0+Count_cf))
    
    # pyomo model instance
    instance = model_construct_inter(data_stru, mode, theta, sim_opt, B_form)
    instance.beta_0.fixed=True
    instance.beta_1.fixed=True
    if B_form>1:
        instance.beta_2.fixed=True

    if sim_opt:
        instance.Obj_1 = Objective(expr = 1)
        instance.Obj = Expression(rule=obj_rule)
    else:
        instance.Obj = Objective(rule=obj_rule, sense=minimize)

    #Try initialize
    try:
        #Simulate the model using scipy
        sim = Simulator(instance, package='casadi') 
        tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
        #Discretize the model using Orthogonal Collocation
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
        #Initialize the discretized model using the simulator profiles
        sim.initialize_model()
    except:
        try:
            #Relax constraint and try intilization again
            instance.ode_mF.deactivate()
            #Simulate the model using scipy
            sim = Simulator(instance, package='casadi') 
            tsim, profiles = sim.simulate(numpoints=300, integrator='idas')
            #Discretize the model using Orthogonal Collocation
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
            #Initialize the discretized model using the simulator profiles
            sim.initialize_model()
            instance.ode_mF.activate()
        except:
            TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')

    solver = SolverFactory('ipopt')
    solver.options["linear_solver"] = "ma97"

    results = solver.solve(instance,tee=True)
    assert results.solver.termination_condition == TerminationCondition.optimal, "Optimization fail"
    results.write()
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    #solver.solve(instance,tee=True).write()
    #instance.display()

    if mode !='DATA':
        fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
    else:
        fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
    sim_inter = inter_model(data_stru, sim_stru, fit_stru)

    if LOUD:
        time = []
        cIn = []
        cH = []
        Jw = []
        Js = []       
        for i in range(data_stru['data_config']['n']):   
            time.extend(sim_stru[i]['time'])
            cIn.extend(sim_stru[i]['cIn'])
            cH.extend(sim_stru[i]['cH'])
            Jw.extend(sim_stru[i]['Jw'])
            Js.extend(sim_stru[i]['Js'])
        
        sim_data = {'time': time,
                'cIn': cIn,
                'cH': cH,
                'Jw': Jw,
                'Js': Js}
        print(sim_data)
        fname = 'sim_data-dat'+str(data_stru['dataset'])       
        #create data frame from dictionary
        sim_datapd = pd.DataFrame(sim_data)
        
        #save dataframe to csv file
        sim_datapd.to_csv(fname+".csv", index=False)
        
        #validate the csv file by importing it
        #print(pd.read_csv(fname+".csv"))

    return fit_stru, sim_stru, sim_inter


# Standard lag diafiltration experiment
mode = 'Lag'
# A.0
theta = {'Lp': 11,
  'beta_c': 15,
  'beta_0': 1,
  'beta_1': 0.01,
  'sigma': 1.0,
  'S0': 0}
data_stru = loadmat('data_library/data_stru-dataset270511.123.mat')['data_stru']
fit_stru_base, sim_stru, sim_inter = solve_model(data_stru, mode, theta, sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

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
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# C.1
data_stru = loadmat('data_library/data_stru-dataset270511.421.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
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
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)

# E.2
data_stru = loadmat('data_library/data_stru-dataset270511.821.mat')['data_stru']
fit_stru, sim_stru, sim_inter = solve_model_B_fix(data_stru, mode, theta=fit_stru_base['parameters'], sim_opt=False, B_form=1)
plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=True)
