%%
% No concentration polarization
% Make folders
if ~exist('fit_B_per_vial\no_concentration_polarization', 'dir')
    mkdir fit_B_per_vial\no_concentration_polarization
end
% mod = {'1cv','1cvmv','101cv','101cvmv','2cv','2cvmv'};
mod = {'1cvmv','101cvmv','2cvmv'};
for i = 1:length(mod)
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{i}];
    if ~exist(dirname, 'dir')
        mkdir(dirname)
    end    
end
%%
% Generate B dependence on retentate concentration plot for multicomponent
% experiments
clear all
Dat = [3.1;3.2];
mod = {'1cv','1cvmv','101cv','101cvmv','2cv','2cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru = load_data(Dat(i));        
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru);
        model_stru = update_model_config(data_stru,model_stru);
        model_stru = update_model_config_pervial(data_stru,model_stru);

        weight_stru = weight_config(best_weight, data_stru, model_stru);
        fit_stru = fitting(reg_method, scale_opt, weight_stru.default, data_stru, model_stru);

        plot_sim(data_stru, model_stru, fit_stru.sim_stru)
        plot_pervial(data_stru, model_stru, fit_stru)
        plot_pervial_vsclassical(data_stru, model_stru, fit_stru)
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname)   
    end
end 

%%
% generate aggeregate B dependence on retentate concentration plots for single component
% experiments
% KCl
clear all
Dat = [1.1;1.2;1.3;1.4;11.1;11.2];
mod = {'1cvmv','101cvmv','2cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru = update_model_config(data_stru(i),model_stru);
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

        plot_sim(data_stru(i), model_stru, fit_stru(i).sim_stru)
        plot_pervial(data_stru(i), model_stru, fit_stru(i))
        plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
    
    plot_aggregate_pervial(data_stru, fit_stru)
    plotname = ['B vs cf_aggregated_',data_stru(1).data_config.namec,'.png'];
    movefile(plotname,dirname)         
end
%%
% MgCl2
clear all
Dat = [2.1;2.2;2.3;12.1;12.2];
mod = {'1cvmv','101cvmv','2cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru = update_model_config(data_stru(i),model_stru);
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

        plot_sim(data_stru(i), model_stru, fit_stru(i).sim_stru)
        plot_pervial(data_stru(i), model_stru, fit_stru(i))
        plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
    
    plot_aggregate_pervial(data_stru, fit_stru)
    plotname = ['B vs cf_aggregated_',data_stru(1).data_config.namec,'.png'];
    movefile(plotname,dirname)         
end



%%
% Considering concentration polarization
% Make folders
if ~exist('fit_B_per_vial\concentration_polarization', 'dir')
    mkdir fit_B_per_vial\concentration_polarization
end
% mod = {'1cv','1cvmv','101cv','101cvmv','2cv','2cvmv'};
mod = {'1cvmv','101cvmv','2cvmv'};
for i = 1:length(mod)
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{i}];
    if ~exist(dirname, 'dir')
        mkdir(dirname)
    end    
end
%%
% Generate B dependence on retentate concentration plot for multicomponent
% experiments
clear all
Dat = [3.1;3.2];
mod = {'1cvmv','101cvmv','2cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru = load_data(Dat(i));        
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru);
        model_stru.concpolar = true;
        model_stru = update_model_config(data_stru,model_stru);
        model_stru = update_model_config_pervial(data_stru,model_stru);

        weight_stru = weight_config(best_weight, data_stru, model_stru);
        fit_stru = fitting(reg_method, scale_opt, weight_stru.default, data_stru, model_stru);

        plot_sim(data_stru, model_stru, fit_stru.sim_stru)
        plot_pervial(data_stru, model_stru, fit_stru)
        plot_pervial_vsclassical(data_stru, model_stru, fit_stru)
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname)   
    end
end 

%%
% generate aggeregate B dependence on retentate concentration plots for single component
% experiments
% KCl
clear all
Dat = [1.1;1.2;1.3;1.4;11.1;11.2];
% mod = {'1cvmv','101cvmv','2cvmv'};
mod = {'101cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru.concpolar = true;
        model_stru = update_model_config(data_stru(i),model_stru);
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

        plot_sim(data_stru(i), model_stru, fit_stru(i).sim_stru)
        plot_pervial(data_stru(i), model_stru, fit_stru(i))
        plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
    
    plot_aggregate_pervial(data_stru, fit_stru)
    plotname = ['B vs cf_aggregated_',data_stru(1).data_config.namec,'.png'];
    movefile(plotname,dirname)         
end
%%
% MgCl2
clear all
Dat = [2.1;2.2;2.3;12.1;12.2];
% mod = {'1cvmv','101cvmv','2cvmv'};
mod = {'101cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru.concpolar = true;
        model_stru = update_model_config(data_stru(i),model_stru);
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

        plot_sim(data_stru(i), model_stru, fit_stru(i).sim_stru)
        plot_pervial(data_stru(i), model_stru, fit_stru(i))
        plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
    
    plot_aggregate_pervial(data_stru, fit_stru)
    plotname = ['B vs cf_aggregated_',data_stru(1).data_config.namec,'.png'];
    movefile(plotname,dirname)         
end

%%
% No concentration polarization
% Make folders
if ~exist('fit_B_per_vial\no_concentration_polarization', 'dir')
    mkdir fit_B_per_vial\no_concentration_polarization
end

mod = {'201cvmv'};
for i = 1:length(mod)
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{i}];
    if ~exist(dirname, 'dir')
        mkdir(dirname)
    end    
end
%%
% generate aggeregate B dependence on retentate concentration plots for single component
% experiments
% KCl
clear all
Dat = [301.1;411.1];
mod = {'201cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\no_concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru = update_model_config(data_stru(i),model_stru);
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

        plot_pervial(data_stru(i), model_stru, fit_stru(i))
        plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
end

%%
% Considering concentration polarization
% Make folders
if ~exist('fit_B_per_vial\concentration_polarization', 'dir')
    mkdir fit_B_per_vial\concentration_polarization
end

mod = {'201cvmv'};
for i = 1:length(mod)
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{i}];
    if ~exist(dirname, 'dir')
        mkdir(dirname)
    end    
end

%%
% generate aggeregate B dependence on retentate concentration plots for single component
% experiments
% KCl
clear all
% Dat = [301.1;501.1;411.1];
Dat = [302.1;402.1];
mod = {'201cvmv'};

reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for j = 1:length(mod)
    legacy_model_name = mod{j};
    dirname = ['fit_B_per_vial\concentration_polarization\model',mod{j}];
    for i = 1:length(Dat)
        data_stru(i) = load_data(Dat(i));  
        % create model structure with default options
        model_stru = model_config(legacy_model_name,data_stru(i));
        model_stru.concpolar = true;
        model_stru = update_model_config(data_stru(i),model_stru);
        
        cF_table = readtable([data_stru(i).filename,'_cF.csv']);
        if cF_table.Properties.VariableNames{1} == "numberOfVial"
            cF = cF_table{:,:};
            model_stru.initialization.C_F0 = cF(1,2:end);
        end
        model_stru = update_model_config_pervial(data_stru(i),model_stru);

        weight_stru = weight_config(best_weight, data_stru(i), model_stru);
        fit_stru(i) = fitting(reg_method, scale_opt, weight_stru.default, data_stru(i), model_stru);

%         plot_pervial(data_stru(i), model_stru, fit_stru(i))
%         plot_pervial_vsclassical(data_stru(i), model_stru, fit_stru(i))
        
        % move folder contains all results to target folder
        movefile([model_stru.filenamestr,'*'],dirname) 
    end 
    plot_aggregate(data_stru, fit_stru)
    plotname = ['B vs t_aggregated_',data_stru(1).data_config.namec,'.png'];
    movefile(plotname,dirname)        
end

function plot_aggregate(data_stru, fit_stru)
% Visualize the aggregate B dependence on retentate concentration from multiple experiments

figure
hold on
% plot B per vial vs retentate concentration
for i = 1:length(fit_stru)
    for j = 1:data_stru(i).data_config.n
        time(j) = fit_stru(i).sim_stru(j).time(end);
        B(j) = fit_stru(i).B{j};
    end    
    h(i) = plot(time, B,...
             's','MarkerSize',10);
    set(h(i), 'markerfacecolor', get(h(i), 'color'))
    
    % Delete the paths, leaving only file names. 
    match = "../data_library/";
    exp_name{i} = erase(data_stru(i).filename,match);   
    exp_name{i} = replace(exp_name{i},'_',' ');
end

xlabel('Time [s]','FontSize',15)
ylabel('B [micrometers/s]','FontSize',15);
title({'\fontsize{15} Salt Permeability vs Retentate Concentration'});

legend(h,exp_name,'Location','eastoutside'); 
set(gcf,'Units','Inches','Position',[0.01 0.01 14 6])
saveas(gcf,['B vs t_aggregated_',data_stru(1).data_config.namec,'.png'])
end