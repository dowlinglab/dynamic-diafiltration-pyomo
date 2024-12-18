import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib import patheffects
import pandas as pd

def loadmat(filename):
    '''
    Read in nested structure(mat file) generated from MATLAB and output dictionaries.

    This function is called instead of using spio.loadmat directly as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries which are still mat-objects.
    Adapted from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    '''
    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mio5_params.mat_struct):
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
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
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
            if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def plot_sim_comparison(data_stru,fit_stru,plot_pred=True,cond=True,lg=False,preface=False):
    '''
    Plot simulation results comparing with measurements
    '''

    t_delay = data_stru['data_raw'][0]['time'][0]

    # plot mass data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    for i in range(data_stru['data_config']['n']):
        plt.plot((data_stru['data_raw'][i]['time']-t_delay)/60,
                 data_stru['data_raw'][i]['mass'],'r.',markersize=4)
        if plot_pred:
            plt.plot((fit_stru['sim_stru'][i]['time']-t_delay)/60,
                     fit_stru['sim_stru'][i]['mV'],'b',linewidth=3,
                     alpha=.6)

    # ghost point for legend
    plt.plot([],[],'r.',markersize=4,label='Measurements')
    if plot_pred:
        plt.plot([],[],'b',linewidth=3,alpha=.6,label='Predictions')

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Mass',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')
    #plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)

    if lg:
        plt.legend(fontsize=10,loc='best')
    plt.show()
    if preface:
        fname = 'mass_preface-dat'+str(data_stru['dataset'])
    else:
        fname = 'mass-dat'+str(data_stru['dataset'])
    fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    plt.plot((data_stru['data_raw'][0]['time'][0]-t_delay)/60,
             fit_stru['sim_stru'][0]['cF'][0],
             'ms',markersize=8,clip_on=False)

    # plot ghost points for legend
    plt.plot([],[],'ms',markersize=8,clip_on=False,label='Retentate (Measurement)')
    if plot_pred:
        plt.plot([],[],'g',linewidth=3,label='Retentate (Prediction)')

    plt.plot([],[],'cs',markersize=8,label='Vial (Measurement)')
    if plot_pred:
        plt.plot([],[],'r^',markersize=8,label='Vial (Prediction)')
        plt.plot([],[],'r-',linewidth=3,alpha=.6,label='Permeate (Prediction)')


    for i in range(data_stru['data_config']['n']):
        plt.plot((data_stru['data_raw'][i]['time'][-1]-t_delay)/60,
                 data_stru['data_raw'][i]['cV_avg'],
                 'cs',markersize=8)
        if cond:
            if isinstance(data_stru['data_raw'][i]['cF_exp'],float):
            	plt.plot((data_stru['data_raw'][i]['time'][-1]-t_delay)/60,
                     	 data_stru['data_raw'][i]['cF_exp'],'ms',markersize=8)
            elif len(data_stru['data_raw'][i]['cF_exp'])>1:
            	plt.plot((data_stru['data_raw'][i]['time']-t_delay)/60,
                     	 data_stru['data_raw'][i]['cF_exp'],'ms',markersize=8)
        if plot_pred:
            plt.plot((fit_stru['sim_stru'][i]['time']-t_delay)/60,
                     fit_stru['sim_stru'][i]['cF'],
                     'g',linewidth=3,alpha=.6)
            plt.plot((fit_stru['sim_stru'][i]['time']-t_delay)/60,
                     fit_stru['sim_stru'][i]['cH'],
                     'r-',linewidth=3,alpha=.6)
            if not preface:
                plt.plot((fit_stru['sim_stru'][i]['time'][-1]-t_delay)/60,
                         fit_stru['sim_stru'][i]['cV'][-1],
                         'r^',markersize=8,alpha=.6)
    if not cond:
        plt.plot((data_stru['data_raw'][-1]['time'][-1]-t_delay)/60,
                 data_stru['data_raw'][-1]['cF_exp'],'ms',markersize=8)

    if preface:
        plt.xlabel('Time',fontsize=24,fontweight='bold')
        plt.ylabel('Concentration',fontsize=24,fontweight='bold')
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
        '''
        plt.text((fit_stru['sim_stru'][7]['time'][-1]-t_delay)/60*0.9,
                 fit_stru['sim_stru'][7]['cF'][-1],
                 'Retentate',fontsize=20,fontweight='bold',
                 horizontalalignment='right',verticalalignment="bottom")

        plt.text((data_stru['data_raw'][8]['time'][-1]-t_delay)/60*0.9,
                 data_stru['data_raw'][8]['cV_avg']*1.3,
                 'Permeate',fontsize=20,fontweight='bold',
                 horizontalalignment='right',verticalalignment="bottom")
        '''
        plt.annotate('Retentate', xy=((fit_stru['sim_stru'][5]['time'][-1]-t_delay)/60,
                     fit_stru['sim_stru'][5]['cF'][-1]),  xycoords='data',
                     xytext=(36, 50), weight='bold', textcoords='offset points',
                     size=20, ha='right', va="center",
                     bbox=dict(boxstyle="round", color = "green", alpha=0.1),
                     arrowprops=dict(arrowstyle="wedge,tail_width=0.3",color = "green",alpha=0.1));

        plt.annotate('Permeate', xy=((data_stru['data_raw'][7]['time'][-1]-t_delay)/60,
                 data_stru['data_raw'][7]['cV_avg']),  xycoords='data',
                     xytext=(36, 50), weight='bold', textcoords='offset points',
                     size=20, ha='right', va="center",
                     bbox=dict(boxstyle="round", color = "red", alpha=0.1),
                     arrowprops=dict(arrowstyle="wedge,tail_width=0.3",color = "red",alpha=0.1));
    else:
        plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
        plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    #plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.minorticks_on()
    plt.tick_params(direction="in",top=True, right=True)
    #plt.tick_params(which="minor",direction="in",top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)

    if lg:
        plt.legend(fontsize=10,loc='best')
    plt.show()
    if preface:
        fname = 'concentration_preface-dat'+str(data_stru['dataset'])
    else:
        fname = 'concentration-dat'+str(data_stru['dataset'])

    fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

def plot_contour(df,show_title=False,preface=False):
    '''
    Plot contours of inividual objectives, two types included:
        1 - fixed sigma,Lp vs. B
        2 - fixed B,Lp vs. sigma
    '''

    f_m = df.Obj_mass.values
    ind_m = np.argmin(f_m)
    f_pc = df.Obj_concentration.values
    ind_pc = np.argmin(f_pc)
    f_rc = df.Obj_retentate_concentration.values
    ind_rc = np.argmin(f_rc)

    F_m = np.reshape(f_m, (50,50))
    F_pc = np.reshape(f_pc, (50,50))
    F_rc = np.reshape(f_rc, (50,50))

    if 'B' in df:
        xx = df.B.values
        if preface:
            xlabelstr = 'B'
            axfontsize = 24
            fixstr = 'preface_fixsig'
        else:
            xlabelstr = 'B [$\mathbf{\mu}$m $\mathbf{\cdot}$ s$\mathbf{^{-1}}$]'
            axfontsize = 16
            fixstr = 'fixsig'
    else:
        xx = df.sigma.values
        if preface:
            xlabelstr = '$\mathbf{\sigma}$'
            axfontsize = 24
            fixstr = 'preface_fixB'
        else:
            xlabelstr = '$\mathbf{\sigma}$ [dimensionless]'
            axfontsize = 16
            fixstr = 'fixB'

    if preface:
        ylabelstr = 'L$\mathbf{_p}$'
    else:
        ylabelstr = r'L$\mathbf{_p}$ [L $\mathbf{ \cdot}$ m$\mathbf{^{-2} \cdot}$h$\mathbf{^{-1} \cdot}$bar$\mathbf{^{-1}}$]'

    yy = df.Lp.values
    X = np.reshape(xx, (50,50))
    Y = np.reshape(yy, (50,50))

    # mass contour
    fig = plt.figure(1,figsize=(4,4))
    cp = plt.contour(X,Y,F_m,10,linewidths=2)
    plt.clabel(cp,cp.levels[::2],inline=True,fontsize=12,colors ='k',fmt='%1.1f')
    plt.plot(xx[ind_m],yy[ind_m],'^',markersize=12,
            markeredgecolor='red',markerfacecolor=[1, .6, .6],clip_on=False)
    masstitle = 'Log$\mathbf{_e}$ transformed \n Mass Objective'
    if show_title:
        plt.title(masstitle,fontsize=16,fontweight='bold')
    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig('contour_'+fixstr+'-mass.png',dpi=300,bbox_inches='tight')

    # permeate concentration contour
    fig = plt.figure(2,figsize=(4,4))
    cp = plt.contour(X,Y,F_pc,10,linewidths=2)
    plt.clabel(cp,cp.levels[::2],inline=True,fontsize=12,colors ='k',fmt='%1.1f')
    plt.plot(xx[ind_pc],yy[ind_pc],'^',markersize=12,
            markeredgecolor='red',markerfacecolor=[1, .6, .6],clip_on=False)
    permtitle = 'Log$\mathbf{_e}$ transformed \n Permeate Concentration Objective'
    if show_title:
       plt.title(permtitle,fontsize=16,fontweight='bold')
    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig('contour_'+fixstr+'-permeate_conc.png',dpi=300,bbox_inches='tight')

    # retentate concentration  contour
    fig = plt.figure(3,figsize=(4,4))
    cp = plt.contour(X,Y,F_rc,10,linewidths=2)
    plt.clabel(cp,cp.levels[::2],inline=True,fontsize=12,colors ='k',fmt='%1.1f')
    plt.plot(xx[ind_rc],yy[ind_rc],'^',markersize=12,
            markeredgecolor='red',markerfacecolor=[1, .6, .6],clip_on=False)
    rettitle = 'Log$\mathbf{_e}$ transformed \n Retentate Concentration Objective'
    if show_title:
       plt.title(rettitle,fontsize=16,fontweight='bold')
    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    if preface:
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.gca().axes.yaxis.set_ticklabels([])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    fig.savefig('contour_'+fixstr+'-retentate_conc.png',dpi=300,bbox_inches='tight')


def plot_sim(sim_stru,color,linetype):
    '''
    Plot simulation results
    '''

    t_delay = sim_stru[0]['time'][0]
    # plot mass prediction
    plt.figure(1,figsize=(4,4))
    for i in sim_stru:
        plt.plot((i['time']-t_delay)/60,i['mV'],color,linestyle=linetype,linewidth=2,alpha=.8)

    # plot retentate concentration prediction
    plt.figure(2,figsize=(4,4))
    for i in sim_stru:
        plt.plot((i['time']-t_delay)/60,i['cF'],color,linestyle=linetype,linewidth=2,alpha=.8)

    # plot permeate concentration prediction comparison
    plt.figure(3,figsize=(4,4))
    for i in sim_stru:
        plt.plot((i['time']-t_delay)/60,i['cH'],color,linestyle=linetype,linewidth=2,alpha=.8)


def plot_sim_show(sig,colors):
    '''
    Figure setup for different sigma values
    '''
    custom_lines = [Line2D([0], [0], color=colors[0],ls='--', lw=3),
                    Line2D([0], [0], color=colors[1],ls='-', lw=3),
                    Line2D([0], [0], color=colors[2],ls=':', lw=3)]
    legendname = ['$\sigma$ = '+str(sig[0]),'$\sigma$ = '+str(sig[1]),'$\sigma$ = '+str(sig[2])]

    fig = plt.figure(1)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Mass [g]',fontsize=16,fontweight='bold')
    #plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    #plt.xlim(left=0)
    plt.ylim(bottom=0)
    fig.savefig('sigma_sensitivity-mass.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(2)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Retentate Conc. [mM]',fontsize=16,fontweight='bold')
    #plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
    plt.legend(custom_lines,legendname,fontsize=15,loc='best')
    fig.savefig('sigma_sensitivity-reten_conc.png',dpi=300,bbox_inches='tight')

    fig = plt.figure(3)
    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate Conc. [mM]',fontsize=16,fontweight='bold')
    #plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
    fig.savefig('sigma_sensitivity-perme_conc.png',dpi=300,bbox_inches='tight')

def plot_contour_sig_sen(data_stru, contour_sig_stru, Vmin, Vmax, Level, Colorbar_ticks, Manual_locations, name_append, filled=False, Jw_filter=True, bar=False):
    '''
    Plot contours for sigma sensitivity (difference in predictions) over different experiment conditions.

    '''
    diaf_mode = False
    if 'cf0' in contour_sig_stru:
        x = contour_sig_stru['cf0']
        xlabelstr = '$\mathbf{c_f}$(t=0) [mM]'
    else:
        x = contour_sig_stru['cd']
        xlabelstr = '$\mathbf{c_d}$ [mM]'
        diaf_mode = True
    y = contour_sig_stru['delp']
    ylabelstr = '$\mathbf{\Delta}$P [psi]'

    X,Y = np.meshgrid(x, y)
    axfontsize = 16

    # gas constant
    R = 8.314e-5 #[cm^3 bar / micromol / K]
    # temperature
    T = data_stru['data_config']['Temp'] #[K]
    # number of dissolved species
    ni = data_stru['data_config']['ni']

    if not diaf_mode:
        # Jw0 = delp - sigma * niRT * cf0
        Jw0 = Y/14.504 - 1*ni*R*T*X
        jw0 = 1*ni*R*T*np.array(x)*14.504
    else:
        cf0 = 5
        Jw0 = Y/14.504 - 1*ni*R*T*cf0
        jw0 = 1*ni*R*T*cf0*14.504

    def patch_back():
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        xy = (xmin,ymin)
        width = xmax - xmin
        height = ymax - ymin

        # create the patch and place it in the back of countourf (zorder!)
        p = patches.Rectangle(xy, width, height, color ='gray', zorder=-10)
        ax.add_patch(p)


    # difference in mass
    fig = plt.figure(1,figsize=(4,4))
    if filled:
        plt.contourf(X,Y,contour_sig_stru['range_mV'],50,cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X,Y,contour_sig_stru['range_mV'],Level[0],linewidths=3,colors ='k')
    else:
        cp = plt.contour(X,Y,contour_sig_stru['range_mV'],Level[0],linewidths=2)
    manual_locations = Manual_locations[0]
    plt.clabel(cp,inline=True,manual=manual_locations,
               fontsize=15,colors ='k',fmt='%1.1f',zorder=2)

    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections,
                 path_effects=[patheffects.withTickedStroke(angle=300, length=2)])

    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom= y[0])
    fig.savefig('contour_sensitivity'+name_append+'-mass.png',dpi=300,bbox_inches='tight')

    # difference in retentate
    fig = plt.figure(2,figsize=(4,4))
    if filled:
        plt.contourf(X,Y,contour_sig_stru['range_cF'],100,cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X,Y,contour_sig_stru['range_cF'],Level[1],linewidths=3,colors ='k')
    else:
        cp = plt.contour(X,Y,contour_sig_stru['range_cF'],Level[1],linewidths=2)
    manual_locations = Manual_locations[1]
    plt.clabel(cp,inline=True,manual=manual_locations,
               fontsize=15,colors ='k',fmt='%1.1f',zorder=2)

    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections,
                 path_effects=[patheffects.withTickedStroke(angle=300, length=2)])

    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom= y[0])
    fig.savefig('contour_sensitivity'+name_append+'-retentate_conc.png',dpi=300,bbox_inches='tight')

    # difference in permeate
    fig = plt.figure(3,figsize=(4,4))
    if filled:
        plt.contourf(X,Y,contour_sig_stru['range_cH'],50,cmap='spring', vmin=Vmin, vmax=Vmax)
        cp = plt.contour(X,Y,contour_sig_stru['range_cH'],Level[2],linewidths=3,colors ='k')
    else:
        cp = plt.contour(X,Y,contour_sig_stru['range_cH'],Level[2],linewidths=2)
    manual_locations = Manual_locations[2]
    plt.clabel(cp,inline=True,manual=manual_locations,
               fontsize=15,colors ='k',fmt='%1.1f',zorder=2)

    if filled:
        if Jw_filter:
            patch_back()
        else:
            plt.fill_between(x, y[0], jw0, color='gray', zorder=2)
    else:
        cg = plt.contour(X, Y, Jw0, [0], colors='orangered')
        plt.setp(cg.collections,
                 path_effects=[patheffects.withTickedStroke(angle=300, length=2)])

    plt.xlabel(xlabelstr,fontsize=axfontsize,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=axfontsize,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(bottom= y[0])
    fig.savefig('contour_sensitivity'+name_append+'-permeate_conc.png',dpi=300,bbox_inches='tight')

    # uniform colorbar
    if filled & bar:
        cmap = mpl.cm.spring
        norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        # horizontal
        fig, ax = plt.subplots(figsize=(12, 1))
        fig.subplots_adjust(bottom=0.5)
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal',
                                        ticks=Colorbar_ticks,
                                        extend='max')
        ax.tick_params(labelsize = 15)
        fig.savefig('colorbar'+name_append+'-horizontal.png',dpi=300,bbox_inches='tight')
        #vertical
        fig, ax = plt.subplots(figsize=(1, 4))
        fig.subplots_adjust(left=0.5)
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical',
                                        ticks=Colorbar_ticks,
                                        extend='max')
        ax.tick_params(labelsize = 15)
        fig.savefig('colorbar'+name_append+'-vertical.png',dpi=300,bbox_inches='tight')

def plot_conc_range(df_f,df_d,lg=True):
    '''
    Plot experiment space - permeate concentration vs. retentate concentration
    '''
    fig = plt.figure(figsize=(4,4))
    # plot filtration
    for i in range(len(df_f.columns)//2):
        xstr = 'F'+str(i+1)+'_cf'
        ystr = 'F'+str(i+1)+'_cp'
        plt.plot(df_f[xstr],df_f[ystr],'^',markersize=8,alpha=.8,clip_on=False)

    # plot diafiltration
    for i in range(len(df_d.columns)//2):
        xstr = 'D'+str(i+1)+'_cf'
        ystr = 'D'+str(i+1)+'_cp'
        plt.plot(df_d[xstr],df_d[ystr],'s',markersize=8,alpha=.8,clip_on=False)

        # ghost point for legend
    plt.plot([],[],'k^',markersize=8,markerfacecolor='white',label='Filtration')
    plt.plot([],[],'ks',markersize=8,markerfacecolor='white',label='Diafiltration')

    plt.xlabel('Retentate [mM]',fontsize=16,fontweight='bold')
    plt.ylabel('Permeate [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    if lg:
        plt.legend(fontsize=15,loc='best')
    plt.show()
    fname = 'concentration_range'
    fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

def plot_cr_measure(data_stru,fit_stru,cr_pred,ybottom,lg=False,cond=True):
    '''
    Plot retentate concentration vs. time from measurements and mass balance
    '''

    t_delay = data_stru['data_raw'][0]['time'][0]

    # plot retentate concentration data
    fig = plt.figure(figsize=(4,4))
    if cond:
        plt.plot((data_stru['data_raw'][0]['time'][0]-t_delay)/60,
                 fit_stru['sim_stru'][0]['cF'][0],
                 'ms',markersize=8,clip_on=False)
    else:
        plt.plot((data_stru['data_raw'][0]['time'][0]-t_delay)/60,
                 data_stru['data_config']['C_F0'],
                 'go',markersize=8,clip_on=False)


    for i in range(data_stru['data_config']['n']):
        if cond:
            if isinstance(data_stru['data_raw'][i]['cF_exp'],float):
                plt.plot((data_stru['data_raw'][i]['time'][-1]-t_delay)/60,
                         data_stru['data_raw'][i]['cF_exp'],'ms',markersize=8)
            elif len(data_stru['data_raw'][i]['cF_exp'])>1:
                plt.plot((data_stru['data_raw'][i]['time']-t_delay)/60,
                         data_stru['data_raw'][i]['cF_exp'],'ms',markersize=8)

    # plot retentate concentration calculated from mass balance comparison
    if len(cr_pred) >0:
        for i in range(data_stru['data_config']['n']):
            plt.plot((data_stru['data_raw'][i]['time'][-1]-t_delay)/60,
                     cr_pred[i],'go',markersize=8)

    plt.xlabel('Time [min]',fontsize=16,fontweight='bold')
    plt.ylabel('Concentration [mM]',fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in",top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=ybottom)

    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)

    plt.plot([],[],'ms',markersize=8,label='Retentate \n (Measurements)')
    plt.plot([],[],'go',markersize=8,label='Retentate \n (Calculated)')
    if lg:
        plt.legend(fontsize=15,loc='best')

    fname = 'cr_measure-dat'+str(data_stru['dataset'])
    fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')

def calib_curve_cond(calib_curve):
    '''
    Plot calibration curve for conductivity probe
    '''
    x = calib_curve.Conductivity.values
    y = calib_curve.Concentration.values

    fig = plt.figure(figsize=(4,4))
    plt.plot(x,y,'bo',markersize=8)

    # calculate the trendline
    z = np.polyfit(x, y, 1)
    l = np.poly1d(z)
    correlation = np.corrcoef(x, y)[0,1]
    r_squared = correlation**2

    plt.plot(x,l(x),'b:',linewidth=3,alpha=.7)

    eqn = 'y=%.3fx+%.2f \n R$\mathbf{^{2}}$=%.4f'%(z[0],z[1],r_squared)
    plt.annotate(eqn, xy=(x[4],y[4]),  xycoords='data',
                 xytext=(-30, 60), weight='bold', textcoords='offset points',
                 size=12, ha='center', va="center",
                 bbox=dict(boxstyle="round", color = "b", alpha=0.1),);

    xlabelstr='Conductivity [$\mathbf{\mu}$S $\mathbf{\cdot}$ cm$\mathbf{^{-1}}$]'
    ylabelstr='Concentration [mM]'
    plt.xlabel(xlabelstr,fontsize=16,fontweight='bold')
    plt.ylabel(ylabelstr,fontsize=16,fontweight='bold')
    plt.xticks(fontsize=15,rotation=45)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tick_params(direction="in",top=True, right=True)

    fig.savefig('calib_curve.png',dpi=300,bbox_inches='tight')

def plot_conc_comparison(data_stru_f,fit_stru_f,data_stru_d,fit_stru_d,plot_pred=True):
    '''
    Plot concentration results comparing with measurements
    '''

    t_delay_f = data_stru_f['data_raw'][0]['time'][0]
    t_delay_d = data_stru_d['data_raw'][0]['time'][0]

    # plot concentration data/prediction comparison
    fig = plt.figure(figsize=(4,4))
    plt.plot(data_stru_f['data_raw'][0]['time'][0]-t_delay_f,
             data_stru_f['cF_ICP'][0],'mv',
             markersize=6,clip_on=False)
    plt.plot(data_stru_d['data_raw'][0]['time'][0]-t_delay_d,
             data_stru_d['cF_ICP'][0],'k^',
             markersize=6,clip_on=False)

    plt.plot(data_stru_f['data_raw'][0]['time'][0]-t_delay_f,
             fit_stru_f['sim_stru'][0]['cF'][0],
             'bs',markersize=6,clip_on=False)
    plt.plot(data_stru_d['data_raw'][0]['time'][0]-t_delay_d,
             fit_stru_d['sim_stru'][0]['cF'][0],
             'rs',markersize=6,clip_on=False)

    # plot ghost points for legend
    plt.plot([],[],'kv',markersize=6,clip_on=False,label='Filtration Retentate ICP')
    plt.plot([],[],'bs',markersize=6,clip_on=False,label='Filtration Retentate conductivity')
    plt.plot([],[],'bo',markersize=6,label='Filtration Vial ICP')

    plt.plot([],[],'k^',markersize=6,clip_on=False,label='Diafiltration Retentate ICP')
    plt.plot([],[],'rs',markersize=6,clip_on=False,label='Diafiltration Retentate conductivity')
    plt.plot([],[],'ro',markersize=6,label='Diafiltration Vial ICP')

    if plot_pred:
        plt.plot([],[],'g-.',linewidth=2,label='Filtration Retentate')
        plt.plot([],[],'g',linewidth=2,label='Diafiltration Retentate')
        plt.plot([],[],'k-.',linewidth=2,alpha=.6,label='Filtration Permeate')
        plt.plot([],[],'k',linewidth=2,alpha=.6,label='Diafiltration Permeate')

    for i in range(data_stru_f['data_config']['n']):
        plt.plot(data_stru_f['data_raw'][i]['time'][-1]-t_delay_f,
                 data_stru_f['data_raw'][i]['cV_avg'],
                 'bo',markersize=6)
        plt.plot(data_stru_f['data_raw'][i]['time'][-1]-t_delay_f,
                 data_stru_f['data_raw'][i]['cF_exp'],
                 'bs',markersize=6)
        if plot_pred:
            plt.plot(fit_stru_f['sim_stru'][i]['time']-t_delay_f,
                     fit_stru_f['sim_stru'][i]['cF'],
                     'g-.',linewidth=2)
            plt.plot(fit_stru_f['sim_stru'][i]['time']-t_delay_f,
                     fit_stru_f['sim_stru'][i]['cH'],
                     'k-.',linewidth=2,alpha=.6)

    for i in range(data_stru_d['data_config']['n']):
        plt.plot(data_stru_d['data_raw'][i]['time'][-1]-t_delay_d,
                 data_stru_d['data_raw'][i]['cV_avg'],
                 'ro',markersize=6)
        plt.plot(data_stru_d['data_raw'][i]['time'][-1]-t_delay_d,
                 data_stru_d['data_raw'][i]['cF_exp'],
                 'rs',markersize=6)
        if plot_pred:
            plt.plot(fit_stru_d['sim_stru'][i]['time']-t_delay_d,
                     fit_stru_d['sim_stru'][i]['cF'],
                     'g',linewidth=2)
            plt.plot(fit_stru_d['sim_stru'][i]['time']-t_delay_d,
                     fit_stru_d['sim_stru'][i]['cH'],
                     'k',linewidth=2,alpha=.6)

    plt.plot(data_stru_f['data_raw'][-1]['time'][-1]-t_delay_f,
             data_stru_f['cF_ICP'][-1],'kv',
             markersize=6,clip_on=False)
    plt.plot(data_stru_d['data_raw'][-1]['time'][-1]-t_delay_d,
             data_stru_d['cF_ICP'][-1],'k^',
             markersize=6,clip_on=False)

    plt.xlabel('Time [s]',fontsize=16)
    plt.ylabel('Concentration [mM]',fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.minorticks_on()
    plt.tick_params(direction="in",top=True, right=True)
    plt.tick_params(which="minor",direction="in",top=True, right=True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(fontsize=12,bbox_to_anchor=(1.1, 1.05),loc='upper left')
    plt.show()
    fname = 'concentration'
    fig.savefig(fname+'.png',dpi=300,bbox_inches='tight')
