"""
Library of functions for diafiltration experiment modeling
Xinhong Liu
University of Notre Dame
"""

import numpy as np
import pandas as pd
import scipy.io as spio
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import idaes
import time
import copy
import os
import json
from sklearn.metrics import r2_score

from pyomo.environ import *
from pyomo.dae import *
import idaes.core.util.scaling as iscale


def loadmat(filename):
    '''
    Read in nested structure(mat file) generated from MATLAB and output dictionaries.

    This function is called instead of using spio.loadmat directly as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries which are still mat-objects.
    Adapted from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    
    Arguments:
        filename: str, location + filename
    
    Returns:
        dictionary contains structured data
    '''
    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _tolist(d[key])

        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            elif isinstance(sub_elem, int):
                elem_list.append(float(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    
    return _check_keys(data)


def plot_sim_comparison(data_stru,sim_stru,stirc_mass=False,plot_pred=True,lg=False,LOUD=False):
    '''
    Plot simulation results comparing with measurements
    
    Arguments:
        data_stru: dict, experimental data dictionary
        sim_stru: dict, simulation results dictionary
        stirc_mass: boolean, if plot mass change in stirred cell
        plot_pred: boolean, if plot model predictions/simualtions
        lg: boolean, if plot legends
        LOUD: boolean, if store the figures
    
    Actions:
        create plots and store (optional)
    '''
    t_delay = data_stru['data_raw'][0]['time'][0]
    vial1 = data_stru['data_config']['n_v0']-1
    t_start = 0

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    for i in range(data_stru['data_config']['n']):
        plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                 data_stru['data_raw'][i]['mass'],'r.',markersize=4)
        if plot_pred and i >= data_stru['data_config']['n_v0']-1:
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['mV'],'b',linewidth=2,
                     alpha=.6)

    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Mass in Vial [g]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in")
    plt.ylim(bottom=0)
    #plt.show()
    
    if LOUD:
        fname = 'figures/mass-dat'+str(data_stru['dataset'])
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))

    for i in range(data_stru['data_config']['n']):
        if type(data_stru['data_raw'][i]['cV_avg']) != list:
            plt.plot((float(data_stru['data_raw'][i]['time'][-1])-t_delay-t_start)/60,
                     data_stru['data_raw'][i]['cV_avg'],
                     'cs',markersize=6)
        else:
            plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                     data_stru['data_raw'][i]['cV_avg'],'cs',markersize=6,clip_on=False)
        plt.plot([(float(t)-t_delay-t_start)/60 for t in data_stru['data_raw'][i]['time']],
                 data_stru['data_raw'][i]['cF_exp'],'ms',markersize=6,clip_on=False)
        if plot_pred:
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cF'],
                     'g',linewidth=2)
            plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                     sim_stru[i]['cH'],
                     'r-',linewidth=2,alpha=.6)
            if type(data_stru['data_raw'][i]['cV_avg']) != list:
                plt.plot((sim_stru[i]['time'][-1]-t_start)/60,
                         sim_stru[i]['cV'][-1],
                         'r^',linewidth=2,alpha=.6)
            else:
                time = data_stru['data_raw'][i]['time']-t_delay
                index = ~np.isnan(data_stru['data_raw'][i]['cV_avg'])
                t = [time[i] for i, val in enumerate(index) if val]
                f = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cV'], fill_value='extrapolate')  
                plt.plot([(ti-t_start)/60 for ti in t],
                         f(t),
                         'r^',linewidth=2,alpha=.6)
                
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(direction="in",top=True, right=True)
    plt.ylim(bottom=0)
    #plt.show()
    
    if LOUD:
        fname = 'figures/concentration-dat'+str(data_stru['dataset'])
        fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
        
    # plot mass of stirred cell
    if stirc_mass:
        fig = plt.figure(figsize=(4,4))
        for i in range(data_stru['data_config']['n']):
            if plot_pred:
                plt.plot([(time-t_start)/60 for time in sim_stru[i]['time']],
                         sim_stru[i]['mF'],'b',linewidth=2,
                         alpha=.6)

        # ghost point for legend
        plt.plot([],[],'r.',markersize=4,label='Mass (Measurements)')
        if plot_pred:
            plt.plot([],[],'b',linewidth=2,alpha=.6,label='Mass (Predictions)')

        plt.plot([],[],'ms',markersize=6,clip_on=False,label='Retentate (Measurements)')
        if plot_pred:
            plt.plot([],[],'g',linewidth=2,label='Retentate (Predictions)')
        if type(data_stru['data_raw'][i]['cV_avg']) != list:
            plt.plot([],[],'cs',markersize=6,label='Vial (Measurements)')
        else:
            plt.plot([],[],'cs',markersize=6,label='Vial (Measurements)')
        if plot_pred:
            plt.plot([],[],'r^',markersize=6,label='Vial (Predictions)')        
            plt.plot([],[],'r-',linewidth=2,alpha=.6,label='Permeate (Predictions)')

        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass in Stirred Cell [g]',fontsize=16,fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tick_params(direction="in")
        ytop = max([max(sim_stru[i]['mF']) for i in range(data_stru['data_config']['n'])])
        plt.ylim(bottom=0,top=ytop+0.5)
        if lg:
            plt.legend(fontsize=12.5,loc='best')#bbox_to_anchor=(1.02, 0.3),borderaxespad=0,ncol=3)
        #plt.show()
        if LOUD:
            fname = 'figures/stirc_mass-dat'+str(data_stru['dataset'])
            fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')


def save_model(m, LO=True, B_form='single', LOUD=False):
    '''
    Save model results
    
    Arguments:
        m: pyomo model instance
        LO: boolean, if mode is lag/overflow
        B_form: float/str, different solute permeability coefficient B formula 
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration              
        LOUD: boolean, if print out parameter results
    
    Returns:
        fit_stru: dict, parameter fit results
        sim_stru: dict, model predictions
    '''
    fit_stru = dict()
    fit_stru['parameters'] = dict()
    fit_stru['sideparameters'] = dict()
    
    fit_stru['parameters']['Lp'] = value(m.Lp) 
    
    if B_form=='single':
        fit_stru['parameters']['B'] = value(m.B)
    elif B_form=='pervial':
        fit_stru['parameters']['B'] = dict()
        for i in m.n_vial:
            fit_stru['parameters']['B'][i] = value(m.B[i])
    elif B_form=='convection':   
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        fit_stru['parameters']['beta_1'] = value(m.beta_1)
    else:
        fit_stru['parameters']['beta_0'] = value(m.beta_0)
        if not isinstance(B_form, str) and B_form != 0:
            fit_stru['parameters']['beta_1'] = value(m.beta_1)
        if not isinstance(B_form, str) and B_form > 1:
            fit_stru['parameters']['beta_2'] = value(m.beta_2)
        if not isinstance(B_form, str) and B_form > 2:
            fit_stru['parameters']['beta_3'] = value(m.beta_3)
            
    fit_stru['parameters']['sigma'] = value(m.sigma)
    if LO:
        fit_stru['parameters']['S0'] = value(m.S0)
        fit_stru['parameters']['S'] = value(m.S)   
    fit_stru['Obj'] = value(m.Obj)
    fit_stru['obj_m'] = value(m.obj_m)
    fit_stru['obj_cv'] = value(m.obj_cv)
    fit_stru['obj_cr'] = value(m.obj_cr)
    fit_stru['Obj_tru'] = value(m.obj_tru)
    fit_stru['obj_cr_tru'] = value(m.obj_cr_tru)
    
        
    if LOUD:
        # print parameters
        print('Lp = ',value(m.Lp),' L / m / m / hr / bar')
        if isinstance(B_form, str):
            if B_form=='single':
                print('B = ',value(m.B),' micrometers / s')
            elif B_form=='pervial':    
                print('B = ',[value(m.B[i]) for i in m.n_vial],' micrometers / s')
            elif B_form=='convection':   
                print('D/l = ',value(m.beta_0))
                print('H = ',value(m.beta_1))
        else:
            if B_form > 2:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-2,'+',value(m.beta_2),'* c_in ^',B_form-1,'+',value(m.beta_3),'* c_in ^',B_form + ']')
            elif B_form > 1:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form-1,'+',value(m.beta_2),'* c_in ^',B_form,+']')
            elif B_form == 0:
                print('B = Jw*[',value(m.beta_0), ']')
            else:
                print('B = Jw*[',value(m.beta_0),'+',value(m.beta_1),'* c_in ^',B_form,']')

        print('sigma = ',value(m.sigma),' dimensionless')
        print('cD = ',value(m.cD),' mM')

        if LO:
            print('S0 = ',value(m.S0),' g / s')
            print('S = ',value(m.S),' g / hr')
        print('Obj = ' ,value(m.Obj))
        print('Obj(m) = ' ,value(m.obj_m))
        print('Obj(cv) = ' ,value(m.obj_cv))
        print('Obj(cr) = ' ,value(m.obj_cr))
        if LO:
            print('Obj(truncate) = ' ,value(m.obj_tru))
            print('Obj(cr_truncate) = ' ,value(m.obj_cr_tru))
            print('count(cr0) = ' ,value(m.count_cr0))
        print('count = ' ,value(m.count))
        print('count(m) = ' ,value(m.count_m))
        print('count(cv) = ' ,value(m.count_cv))
        print('count(cr) = ' ,value(m.count_cr))
        print('likelihood1 = ' ,value(m.llh1))
        print('likelihood2 = ' ,value(m.llh2)) 
    
    sim_stru = dict()
    
    for n_vial in m.n_vial:
        time = [t * (m.tf[n_vial]-m.ti[n_vial]) + m.ti[n_vial] for t in m.tau]
        cF = [value(m.cF[n_vial,t]) for t in m.tau]
        cIn = [value(m.cIn[n_vial,t]) for t in m.tau]
        cH = [value(m.cH[n_vial,t]) for t in m.tau]
        mV = [value(m.mV[n_vial,t]) for t in m.tau]
        cV = [value(m.cV[n_vial,t]) for t in m.tau]
        Jw = [value(m.Jw[n_vial,t]) for t in m.tau]
        Js = [value(m.Js[n_vial,t]) for t in m.tau]
        
        sim_stru[n_vial-1] = {'time': time,
                'cF': cF,
                'cIn': cIn,
                'cH': cH,
                'mV': mV,
                'cV': cV,
                'Jw': Jw,
                'Js': Js}
        if LO:
            mF = [value(m.mF[n_vial,t]) for t in m.tau]
            sim_stru[n_vial-1]['mF'] = mF
        
        if isinstance(B_form, str) and 'convection' in B_form:
            B = None       
        elif B_form=='single':
            B = [value(m.B) for t in m.tau]
        elif B_form=='pervial':
            B = [value(m.B[n_vial]) for t in m.tau]
        else:
            B = [value(m.B[n_vial,t]) for t in m.tau]        
        sim_stru[n_vial-1]['B'] = B

    return fit_stru, sim_stru


def inter_model(data_stru, sim_stru, fit_stru):
    '''
    Interpolate model predictions at same time stamps of experimental measurements
    
    Arguments:
        data_stru: dict, experimental data dictionary
        sim_stru: dict, model predictions
        fit_stru: dict, parameter fit results
    
    Returns:
        sim_inter: dict, model predictions for experimental measurements 
    '''
    sim_inter = dict()
    for i in range(data_stru['data_config']['n']):
        t_delay = data_stru['data_raw'][0]['time'][0]
        t_meas = data_stru['data_raw'][i]['time'] - t_delay
        f_m = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['mV'], kind='linear', fill_value='extrapolate')
        if i >= data_stru['data_config']['n_v0']-1:            
            ind_m = ~np.isnan(data_stru['data_raw'][i]['mass'])
            t_m = [t_meas[i] for i, val in enumerate(ind_m) if val]
        else:
            t_m=[]
        if type(data_stru['data_raw'][i]['cV_avg']) == list:
            f_cv = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cV'], kind='linear', fill_value='extrapolate')
            ind_cv = ~np.isnan(data_stru['data_raw'][i]['cV_avg'])
            t_cv = [t_meas[i] for i, val in enumerate(ind_cv) if val]
            cV = f_cv(t_cv)
        else:
            cV = sim_stru[i]['cV'][-1]
        f_cf = interpolate.interp1d(sim_stru[i]['time'],sim_stru[i]['cF'], kind='linear', fill_value='extrapolate')
        ind_cf = ~np.isnan(data_stru['data_raw'][i]['cF_exp'])
        t_cf = [t_meas[i] for i, val in enumerate(ind_cf) if val]
        sim_inter[i] = {'time': t_meas,
                        'mV': f_m(t_m),
                        'cF': f_cf(t_cf),
                        'cV': cV}       

    return sim_inter


def model_construct_inter(data_stru, mode, theta=None, sim_opt=False, B_form='single'):
    '''
    Build diafiltration model in pyomo
    
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
    
    Returns:
        m: pyomo model instance
    '''
    # known parameters
    # applied pressure
    delP = data_stru['data_config']['delP'] #[bar]
    # gas constant
    R = 8.314e-5 #[cm^3 bar / micromol / K]
    # temperature
    T = data_stru['data_config']['Temp'] #[K]
    # membrane area
    Am = data_stru['data_config']['Am'] #[cm^2]
    # density
    rho = data_stru['data_config']['rho'] #[g/cm^3]
    # number of dissolved species
    ni = data_stru['data_config']['ni']
    # kinematic viscosity
    nu = 8.927e-3 #[cm^2/s]
    # diameter of stirred cell
    b = 2.2860 #[cm]
    # average velocity within the system
    v = 350/60 * np.pi * b #[cm/s]
    # diffusion coefficient
    if isinstance(data_stru['data_config']['namec'], str) and 'K' in data_stru['data_config']['namec']:
        D = 1.960e-5 #[cm^2/s]  - K+
    else:
        raise NotImplementedError
    # mass transfer coefficien
    k = 0.23 * v**0.57 * D**0.67 / (nu**0.24 * b**0.43)

    # known constants
    t_delay = data_stru['data_raw'][0]['time'][0]
    M_F0 = data_stru['data_config']['M_F0']
    M_O = data_stru['data_config']['M_O']
    cD = data_stru['data_config']['C_D']
    mH = 0.25 #ml
    C_F0 = data_stru['data_config']['C_F0']
    C_V0 = 1e-6
    
    N_VIAL = data_stru['data_config']['n']
    N_V0 = data_stru['data_config']['n_v0']
    if mode !='DATA':
        N_extra = data_stru['data_config']['n_extra']
        N_H = data_stru['data_config']['n_h']
        N_A = data_stru['data_config']['n_A']       
        if mode == 'Overflow':
            N_A0 = 1
            
    Tauf = 1 #scaled ending time

    # create a model object
    m = ConcreteModel()

    # define the independent variable
    m.n_vial = RangeSet(N_VIAL) # number of vials
    m.tau = ContinuousSet(bounds=(0, Tauf))#scaled time
    
    TF_list = [data_stru['data_raw'][i]['time'][-1]-t_delay for i in range(N_VIAL)]
    TF_dict = dict(zip(m.n_vial,TF_list)) # unscaled time elapse for each vial 
    m.tf = Set(initialize=TF_list)
    
    TI_list = [data_stru['data_raw'][i]['time'][0]-t_delay for i in range(N_VIAL)]
    TI_dict = dict(zip(m.n_vial,TI_list)) # unscaled initial time for each vial
    m.ti = Set(initialize=TI_list)

    # parameter initialization
    param_in = dict()
    if theta is None:
        param_in['Lp'] = data_stru['data_config']['Lp0']
        param_in['B'] = data_stru['data_config']['B0']
        param_in['sigma'] = data_stru['data_config']['sigma0']
        param_in['beta_0'] = param_in['B']*36000/param_in['Lp']/delP #param_in['B']
        param_in['beta_1'] = 1
        param_in['beta_2'] = 0
        param_in['beta_3'] = 0
        if mode == 'Lag':
            param_in['S0'] = 0      
        elif mode == 'Overflow':
            param_in['S0'] = -M_O/10 *3600
    else:
        param_in = theta
    # define model parameters
    if sim_opt:        
        m.cD = cD
        m.C_H0 = 1e-6
        #simulation parameters
        m.Lp = Param(initialize=param_in['Lp'],mutable=True)
        
        if B_form=='single':
            m.B = Param(initialize=param_in['B'],mutable=True)
        elif B_form=='pervial':    
            m.B = Param(m.n_vial,initialize=param_in['B'],mutable=True)
        elif B_form=='convection':    
            m.beta_0 = Param(initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)    
        else:
            m.beta_0 = Param(initialize=param_in['beta_0'],mutable=True)
            if not isinstance(B_form, str):
                if B_form != 0:
                    m.beta_1 = Param(initialize=param_in['beta_1'],mutable=True)
                if B_form > 1:
                    m.beta_2 = Param(initialize=param_in['beta_2'],mutable=True)
                if B_form > 2:
                    m.beta_3 = Param(initialize=param_in['beta_3'],mutable=True) 
                m.B = Var(m.n_vial, m.tau, bounds=(1e-6,50))
        m.sigma = Param(initialize=param_in['sigma'],mutable=True)

    else: 
        #parameter estimation
        m.C_H0 = 1e-6
        m.cD = cD
        m.Lp = Var(bounds=(0.5,50),initialize=param_in['Lp'])

        if B_form=='single':
            m.B = Var(bounds=(1e-6,30),initialize=param_in['B'])
        elif B_form=='pervial':    
            m.B = Var(m.n_vial, bounds=(1e-6,30),initialize=param_in['B'])
        elif isinstance(B_form, str) and 'convection' in B_form:
            m.beta_0 = Var(bounds=(1+1e-6,50),initialize=param_in['beta_0'])
            m.H = Var(m.n_vial, m.tau,initialize=0.5)  
            m.beta_1 = Var(bounds=(0.0,1.0),initialize=0.5)         
        else:    
            m.beta_0 = Var(initialize=param_in['beta_0'])
            m.B = Var(m.n_vial, m.tau)
            if B_form != 0:
                m.beta_1 = Var(bounds=(-20,20),initialize=param_in['beta_1'])
            if B_form > 1:
                m.beta_2 = Var(bounds=(-20,20),initialize=param_in['beta_2'])
            if B_form > 2:
                m.beta_3 = Var(bounds=(-20,20),initialize=param_in['beta_3'])
                
        m.sigma = Var(bounds=(0.0,1.0), initialize=param_in['sigma'])
    
    if mode == 'Lag':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True) 
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S0 = Param(initialize=param_in['S0'])
            # g/hr
            m.S = Var(initialize=M_O/sum(m.tf-m.ti)*3600)
    elif mode == 'Overflow':
        if sim_opt:
            m.S0 = Param(initialize=param_in['S0'],mutable=True)
            m.S = Param(initialize=param_in['S'],mutable=True)
        else:
            # g/s
            m.S = Var(initialize=0)
            # g/hr
            m.S0 = Var(initialize=-M_O/10)
            
    # define dependent variables 
    if mode !='DATA':
        m.mF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=M_F0)
    m.cF = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cIn = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=C_F0)
    m.cH = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.mV = Var(m.n_vial, m.tau, domain=NonNegativeReals,initialize=1e-6)
    m.cVmV = Var(m.n_vial, m.tau,initialize=C_V0 * 1e-6)
    m.cV = Var(m.n_vial,  m.tau, domain=NonNegativeReals,initialize=1e-6)

    # define intermediate variables
    m.Jw = Var(m.n_vial, m.tau)
    m.Js = Var(m.n_vial, m.tau)
    if isinstance(B_form, str) and 'convection' in B_form:
        m.Js_exp = Var(m.n_vial, m.tau, bounds=(1+1e-6,1e4))
    elif B_form=='K':   
        m.Js_exp = Var(m.n_vial, m.tau)

    # define derivatives
    if mode !='DATA':
        m.dmF = DerivativeVar(m.mF,wrt=m.tau)
    m.dcF = DerivativeVar(m.cF,wrt=m.tau)
    m.dcH = DerivativeVar(m.cH,wrt=m.tau)
    m.dmV = DerivativeVar(m.mV,wrt=m.tau)
    m.dcVmV = DerivativeVar(m.cVmV,wrt=m.tau)
    
    # define the differential equation as constraints
    def ode_mF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                # (dmF_dt = 0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            elif N_A0-1 < n <= N_A:
                # (dmF_dt = -S0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (- m.S0  - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
            else:
                # (dmF_dt = S) * tf
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # (dmF_dt = -S0 - Jw * Am * rho) * tf
                return m.dmF[n,t] == (- m.S0 - Am * rho * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf   
            else:
                # (dmF_dt = S) * tf
                return m.dmF[n,t] == m.S / 3600  * (TF_dict[n]-TI_dict[n])/Tauf
    if mode !='DATA':
        m.ode_mF = Constraint(m.n_vial, m.tau, rule=ode_mF_rule)
     
    def ode_cF_rule(m, n, t):
        if mode == 'Overflow':
            if n < N_A0:
                return m.dcF[n,t] == 1 / m.mF[n,t] * (Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf        
            elif N_A0-1 < n <= N_A:  
                # (dcF/dt = (cF - cD) * S0 / mF + Am * rho / mF * (cF * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf                     
            else:
                # (dcF/dt = (cD - cF) * S / mF + Am * rho / mF * (cD * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'Lag':
            if n <= N_A:
                # (dcF/dt = (cF - cD) * S0 / mF + Am * rho / mF * (cF * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cF[n,t] - m.cD ) * m.S0 + Am * rho * (m.cF[n,t] * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf                     
            else:
                # (dcF/dt = (cD - cF) * S / mF + Am * rho / mF * (cD * Jw - Js)) * tf
                return m.dcF[n,t] == 1 / m.mF[n,t] * ((m.cD - m.cF[n,t]) * m.S / 3600 + Am * rho * (m.cD * m.Jw[n,t] - m.Js[n,t])) * (TF_dict[n]-TI_dict[n])/Tauf
        elif mode == 'DATA':
            return m.dcF[n,t] == Am * rho / M_F0 * (m.cD * m.Jw[n,t] - m.Js[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf

    m.ode_cF = Constraint(m.n_vial, m.tau, rule=ode_cF_rule)

    # (dcH_dt = Am * rho / mH *(Js - Jw * cH)) * tf
    def ode_cH_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcH[n,t] == Am * rho / m.mV[n,t] * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
        else:
            return m.dcH[n,t] == Am * rho / mH * (m.Js[n,t] - m.cH[n,t] * m.Jw[n,t]) * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cH = Constraint(m.n_vial, m.tau, rule=ode_cH_rule)
    
    # (dmV_dt = Jw * Am * rho) * tf
    def ode_mV_rule(m, n, t):
        return m.dmV[n,t] == m.Jw[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_mV = Constraint(m.n_vial, m.tau, rule=ode_mV_rule)
    
    # (dmVcV_dt = Jw * Am * rho * cH) * tf
    def ode_cVmV_rule(m, n, t):
        if mode !='DATA' and n <= N_H:
            return m.dcVmV[n,t] == Am * rho * m.Js[n,t] * (TF_dict[n]-TI_dict[n])/Tauf
        else:                
            return m.dcVmV[n,t] == m.Jw[n,t] * m.cH[n,t] * Am * rho * (TF_dict[n]-TI_dict[n])/Tauf
    m.ode_cVmV = Constraint(m.n_vial, m.tau, rule=ode_cVmV_rule)
    
    # intermediate variables constraints
    # mF(end) - mF(0) = m_o
    def eqn_S_rule(m):
        return m.mF[m.n_vial.last(),Tauf] - M_F0 == M_O
    if mode !='DATA' and not sim_opt:
        m.eqn_S = Constraint(rule=eqn_S_rule)    
    # cIn = (cF - cH) .* exp(Jw./ k) + cH
    def eqn_cIn_rule(m, n, t):
        return m.cIn[n,t] == (m.cF[n,t] - m.cH[n,t]) * exp(m.Jw[n,t]/ k) + m.cH[n,t]
    m.eqn_cIn = Constraint(m.n_vial, m.tau, rule=eqn_cIn_rule)
    
    # Jw = Lp*(delP - (cIn - cH).' *(ni.*sigma)*R*T)
    # cm/s
    def eqn_Jw_rule(m, n, t):
        return m.Jw[n,t]*36000 == m.Lp *(delP - (m.cIn[n,t] - m.cH[n,t])*ni*m.sigma*R*T)
    m.eqn_Jw = Constraint(m.n_vial, m.tau, rule=eqn_Jw_rule)
    
    # Js = B*(cIn - cH)
    # micromole/cm^2/s
    def eqn_Js_rule(m, n, t):
        if B_form=='single':
            return m.Js[n,t]*10000 == m.B * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='pervial':    
            return m.Js[n,t]*10000 == m.B[n] * (m.cIn[n,t] - m.cH[n,t])
        elif B_form=='convection':
            return m.Js[n,t] == m.Jw[n,t] * m.H[n,t] * (m.cIn[n,t] * m.Js_exp[n,t] - m.cH[n,t]) / (m.Js_exp[n,t]-1)
        else:
            #return m.Js[n,t]*10000 ==  m.B[n,t] * (m.cIn[n,t] - m.cH[n,t]) # No Jw in Js
            return m.Js[n,t]*10000 ==  (m.Jw[n,t] * m.B[n,t] * 10000) * (m.cIn[n,t] - m.cH[n,t]) #Jw in Js
            
    m.eqn_Js = Constraint(m.n_vial, m.tau, rule=eqn_Js_rule)
    
    def eqn_Js_exp_rule(m, n, t):
        if B_form=='convection':
            return m.Js_exp[n,t]  == exp(m.Jw[n,t]/m.beta_0*10000)   
        else:
            return Constraint.Skip
    m.eqn_Js_exp = Constraint(m.n_vial, m.tau, rule=eqn_Js_exp_rule)
    
    def eqn_H_rule(m,n,t):
        if B_form=='convection':
            return m.H[n,t]  == m.beta_1
        else:
            return Constraint.Skip
    m.eqn_H = Constraint(m.n_vial, m.tau, rule=eqn_H_rule)
    
    # cV = (mV*cV) /mV
    def eqn_cV_rule(m, n, t):
        return m.mV[n,t] * m.cV[n,t] == m.cVmV[n,t]
    m.eqn_cV = Constraint(m.n_vial, m.tau, rule=eqn_cV_rule)# numercally better without division
    
    # link the variables from different vials 
    def mF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.mF[n,Tauf]==m.mF[n+1,0]  
    if mode !='DATA':
        m.mF_linking = Constraint(m.n_vial,rule=mF_linking_rule) 
    
    def cF_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cF[n,Tauf]==m.cF[n+1,0]                                                                  
    m.cF_linking = Constraint(m.n_vial,rule=cF_linking_rule)   

    def cH_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        else:
            return m.cH[n,Tauf]==m.cH[n+1,0]                                                                  
    m.cH_linking = Constraint(m.n_vial,rule=cH_linking_rule)      

    def mV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.mV[n,Tauf]==m.mV[n+1,0]
        else:
            return 1e-6==m.mV[n+1,0]                                                                  
    m.mV_linking = Constraint(m.n_vial,rule=mV_linking_rule) 
    
    def cVmV_linking_rule(m,n):
        if n == m.n_vial.last():
            return Constraint.Skip
        elif mode !='DATA'and N_H < N_V0 and n < N_H:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_H < N_V0 and N_H < n < N_V0-1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        elif mode !='DATA'and N_V0-1 < N_extra and N_V0 <= n < N_extra+1:
            return m.cVmV[n,Tauf]==m.cVmV[n+1,0]
        else:
            return m.cV[n,Tauf]*1e-6==m.cVmV[n+1,0]                                                                  
    m.cVmV_linking = Constraint(m.n_vial,rule=cVmV_linking_rule)  
    
   # B per vial, skip startup vials
    def B_form_rule1(m,n):
        if n < N_V0-1:
            return m.B[n] == m.B[n+1]
        else:
            return Constraint.Skip

    # B-cf dependence
    def B_form_rule2(m,n,t):
        if not isinstance(B_form, str):
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**(B_form - 2) + m.beta_2 * m.cF[n,t]**(B_form - 1) + m.beta_3 * m.cF[n,t]**B_form #+ m.beta_4 * m.cH[n,t]        
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**(B_form - 1) + m.beta_2 * m.cF[n,t]**B_form #+ m.beta_4 * m.cH[n,t]
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cF[n,t]**B_form
        else:
            return Constraint.Skip

    # B-cin dependence
    def B_form_rule3(m,n,t):
        if not isinstance(B_form, str):
            if B_form > 2:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**(B_form - 2) + m.beta_2 * m.cIn[n,t]**(B_form - 1) + m.beta_3 * m.cIn[n,t]**B_form# + m.beta_4 * m.cH[n,t]        
            elif B_form > 1:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**(B_form - 1) + m.beta_2 * m.cIn[n,t]**B_form# + m.beta_4 * m.cH[n,t]
            elif B_form == 0:
                return m.B[n,t] == m.beta_0
            elif B_form < 0:
                return m.B[n,t] * m.cIn[n,t]**(-B_form) == m.beta_0 * m.cIn[n,t]**(-B_form) + m.beta_1
            else:
                return m.B[n,t] == m.beta_0 + m.beta_1 * m.cIn[n,t]**B_form
        else:
            return Constraint.Skip
    
    if B_form=='pervial'and not sim_opt: 
        m.b_form = Constraint(m.n_vial,rule=B_form_rule1)
        pass
    else:
        m.b_form = Constraint(m.n_vial,m.tau,rule=B_form_rule3) 

    #initial conditions
    def _init(m): 
        def firstNonNan(listfloats):
            for item in listfloats:
                if np.isnan(item) == False:
                    return item
        if mode !='DATA':
            yield m.mF[1,0]==data_stru['data_config']['M_F0']
            yield m.cH[1,0]==m.C_H0 #1e-6
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
        else:
            if type(data_stru['data_raw'][0]['cV_avg']) == list:
                C_H0 = firstNonNan(data_stru['data_raw'][0]['cV_avg'])
            else:
                C_H0 = data_stru['data_raw'][0]['cV_avg']
            yield m.cH[1,0]==C_H0 * 0.8
            yield m.cF[1,0]==firstNonNan(data_stru['data_raw'][0]['cF_exp'])
               
        yield m.cVmV[1,0]==1e-6*1e-6
        yield m.mV[1,0]==1e-6
        yield ConstraintList.End
        
    m.con_boundary = ConstraintList(rule=_init)
        
    return m


def solve_model(data_stru, mode, theta=None, sim_opt=False, B_form='single', LOUD=False):
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
    #instance.pprint()
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
        #Discretize the model using finite_difference
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')
        #Initialize the discretized model using the simulator profiles
        sim.initialize_model()
    except Exception as e:
        print(f"Initialization failed: {e}. Applying discretization without initialization.")
        TransformationFactory('dae.finite_difference').apply_to(instance, nfe=300, scheme='BACKWARD')

    solver = SolverFactory('ipopt')
    solver.options["linear_solver"] = "ma97"
    solver.options["max_iter"] = 3000
    #solver.options["halt_on_ampl_error"] = "yes" # option 1
    #solver.options["print_level"] = 1 # option 2
    try:
        results = solver.solve(instance,tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            results.write()
        else:
            print(f"Solver termination: {results.solver.termination_condition}, resolving...")
            results = solver.solve(instance,tee=True)
            assert results.solver.termination_condition == TerminationCondition.optimal, (
                    f"Solver failed: Non-optimal termination condition "
                    f"{results.solver.termination_condition}. ")
    except:
        print("Retrying with adjusted solver settings...")
        solver.options["linear_solver"] = "ma57"
        solver.options['max_iter']=5000       
        results = solver.solve(instance,tee=True)
        assert results.solver.termination_condition == TerminationCondition.optimal, (
                f"Solver failed again: Non-optimal termination condition "
                f"{results.solver.termination_condition}. "
                "Try alternative initialization or further debugging.")

    # Process results if optimal
    if results.solver.termination_condition == TerminationCondition.optimal:
        if mode !='DATA':
            fit_stru, sim_stru=save_model(instance, LO=True, B_form=B_form, LOUD=True)
        else:
            fit_stru, sim_stru=save_model(instance, LO=False, B_form=B_form, LOUD=True)
        sim_inter = inter_model(data_stru, sim_stru, fit_stru)
    else:
        print("Non-optimal solution after retries. Check initialization or parameter settings.")
        fit_stru, sim_stru, sim_inter = None, None, None
        return fit_stru, sim_stru, sim_inter

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


def nested_dict_values(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all values of nested dictionaries
    '''
    dict_value=[]
    # Iterate over all values of given dictionary
    for value in dict_nest.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for v in  nested_dict_values(value):
                dict_value.append(v)
        else:
            # If value is not dict type then append the value
            dict_value.append(value)
    return dict_value

def nested_dict_keys(dict_nest):
    ''' 
    This function accepts a nested dictionary as argument
        and iterates over all keys of nested dictionaries
    '''   
    dict_key=[]
    # Iterate over all values of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):           
            # If value is dict then iterate over all its keys
            for v in nested_dict_keys(dict_nest[key]):
                dict_key.append(key+'-'+str(v))
        else:
            # If value is not dict type then yield the value
            dict_key.append(key)
    return dict_key

def nested_dict_update(dict_nest,new_value,i=0):
    ''' 
    This function accepts a nested dictionary and a new list of values as arguments
        and updates all values of nested dictionaries
    '''    
    # Iterate over all keys of given dictionary
    for key in dict_nest.keys():
        # Check if value is of dict type
        if isinstance(dict_nest[key], dict):
            # If value is dict then iterate over all its keys
            for v in nested_dict_update(dict_nest[key],new_value,i=i):
                v = new_value[i]
                i += 1
        else:
            # If value is not dict type then update the value
            dict_nest[key] = new_value[i]
            i += 1
    return dict_nest


def calc_FIM(data_stru, mode, theta=None, step=1e-8, formula='backward', B_form='single'):
    """
    Calculate FIM
    
    Arguments:
        data_stru: dict, experimental data dictionary
        mode: str, experiment mode, {DATA, lag, overflow}
        theta: dict, preset parameter values
        step: float, relative step change for parameters
        formula: str, finite difference scheme option, {backward, forward, central}
        B_form: float/str, different solute permeability coefficient formula
                'single' - constant B
                'per vial' - discrete B per vial
                'convection' - convection-diffussion model
                0, -0.5, 0.5, 1, 2, 3 - order of dipendence on concentration 
    
    Returns:
        doe_stru: dict, FIM results
    """
    sim_opt = True
    if theta is None:
        sim_opt = False
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta, sim_opt, B_form)
        theta = fit_stru_p['parameters']
        sim_opt = True
    theta_p = copy.deepcopy(theta)
    theta_p_v=nested_dict_values(theta_p)
    theta_p1_v = []
    theta_p2_v = []

    # Step changes
    for i in theta_p_v: 
        
        if formula == 'central':
            if i!=0:
                theta_p1_v.append(i * (1-step))
                theta_p2_v.append(i * (1+step))
            else:
                theta_p1_v.append(-step)
                theta_p2_v.append(step)
        elif formula == 'backward':
            if i!=0:
                theta_p1_v.append(i * (1-step))
            else:
                theta_p1_v.append(-step)
            theta_p2_v.append(i)
        else: #forward
            theta_p1_v.append(i)
            if i!=0:
                theta_p2_v.append(i * (1+step))
            else:
                theta_p2_v.append(step)
    
    # Prepare prediction covariance matrix        
    if sim_opt == True:
        fit_stru_p, sim_stru_p, sim_inter_p = solve_model(data_stru, mode, theta_p, sim_opt, B_form)
    var_pred=[]
    for n_vial in range(data_stru['data_config']['n']):
        var_pred = np.append(var_pred,0.01 ** 2 * np.ones(len(sim_inter_p[n_vial]['mV'])))
        var_pred = np.append(var_pred,(0.03 * sim_inter_p[n_vial]['cV'])**2)
        var_pred = np.append(var_pred,(0.003 * sim_inter_p[n_vial]['cF'])**2)
    cov_pred = np.diag(var_pred)
    
    # Prepare Jacobian
    Jac = []
    for i in range(len(theta_p_v)):
        theta_pk1_v = theta_p_v.copy()
        theta_pk2_v = theta_p_v.copy()
        theta_pk1_v[i] = theta_p1_v[i]    
        theta_pk2_v[i] = theta_p2_v[i]       
        theta_pk1 = copy.deepcopy(theta_p)
        theta_pk2 = copy.deepcopy(theta_p)
        theta_pk1 = nested_dict_update(theta_pk1,theta_pk1_v)
        theta_pk2 = nested_dict_update(theta_pk2,theta_pk2_v)
        print(theta_pk1)
        print(theta_pk2)
        
        sim_opt = True
        fit_stru_pk1, sim_stru_pk1, sim_inter_pk1 = solve_model(data_stru, mode, theta_pk1, sim_opt, B_form)
        fit_stru_pk2, sim_stru_pk2, sim_inter_pk2 = solve_model(data_stru, mode, theta_pk2, sim_opt, B_form)
        jac=[]
        for n_vial in range(data_stru['data_config']['n']):
            jac = np.append(jac,sim_inter_pk2[n_vial]['mV']-sim_inter_pk1[n_vial]['mV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cV']-sim_inter_pk1[n_vial]['cV'])
            jac = np.append(jac,sim_inter_pk2[n_vial]['cF']-sim_inter_pk1[n_vial]['cF'])
        # Predictions in row
        if theta_p_v[i]!=0:
            delta = step*theta_p_v[i]
        else:
            delta = step
        if len(Jac) == 0:
            Jac = jac/delta
        else:
            Jac = np.vstack([Jac, jac/delta])

    doe_stru = dict()
    if formula == 'central':
        doe_stru['Jac'] = Jac/2
    else:
        doe_stru['Jac'] = Jac
    FIM = doe_stru['Jac'] @ np.linalg.inv(cov_pred) @ doe_stru['Jac'].T
    doe_stru['Jac'] = doe_stru['Jac'].tolist()
    doe_stru['FIM'] = FIM.tolist()

    # Compute eigenvalues of FIM
    w, v = np.linalg.eigh(FIM)
    doe_stru['eig_val'] = w.tolist()
    doe_stru['eig_vec'] = v.tolist() # in columns
    theta_p_name=nested_dict_keys(theta_p)
    maxind = np.argmax(abs(v), axis=0)
    doe_stru['eig_dir']=[theta_p_name[i] for i in (maxind)]
    doe_stru['trace'] = np.trace(FIM).tolist()
    doe_stru['det'] = np.linalg.det(FIM).tolist()
    doe_stru['min_eig'] = min(w)
    doe_stru['cond'] = max(w) / min(w)
    try:
        doe_stru['V'] = np.linalg.inv(FIM).tolist()
        doe_stru['std'] = np.sqrt(np.diag(doe_stru['V'])).tolist()
    except:
        print('No inverse of FIM')

    print(doe_stru)
    return doe_stru

def correlation_from_covariance(covariance):
    ''' 
    This function calculate correlation matrix from covariance matrix.
    '''    
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    print(correlation)
    return correlation

def store_json(file_name,structure):
    with open(file_name, "w") as json_file:
        json.dump(structure, json_file)