% Generate heatmap for doe on different experiment conditions
clear all
% Base experiment
Dat = 501.1;
% Experiment conditions
cf0 = 5:1:80;
delp = 10:1:120;
% Parameters evaluated at
theta = [3.8994;0.29489;1];% NF90.5 diafiltration

view_heatmap = true;
Jw_filter = true;

% perform default regression
reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.equal_statistical = true; % Hybrid weight
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

step = 1e-10;
formula = 'forward';

% Run simulation
data_stru = load(['data/data_stru-dataset', num2str(Dat), '.mat']).data_stru; 
legacy_model_name = '201cvmv';

% Initialization
trace = zeros(length(delp),length(cf0));
determinant = zeros(length(delp),length(cf0));
min_eig = zeros(length(delp),length(cf0));
cond = zeros(length(delp),length(cf0));

for p = 1:length(delp)
    delP = delp(p);
    data_stru.data_config.delP = delP/14.504;
    
    model_stru = model_config(legacy_model_name,data_stru);
    model_stru.concpolar = false;
    model_stru = update_model_config(data_stru,model_stru);
    
    for c = 1:length(cf0)
        C_F0 = cf0(c);
        model_stru.initialization.C_F0 = C_F0;
            
        weight_stru = weight_config(best_weight, data_stru, model_stru);
        doe_stru = calc_FIM(scale_opt, weight_stru.equal_statistical, theta, step, formula, data_stru, model_stru);
        trace(p,c) = doe_stru.trace;
        determinant(p,c) = doe_stru.det;
        min_eig(p,c) = doe_stru.min_eig;
        cond(p,c) = doe_stru.cond;
        
        % Jw
        if Jw_filter
            nRT = data_stru.data_config.ni * 8.314e-5 * data_stru.data_config.Temp;
            sim_stru = sim_model(theta, data_stru, model_stru);
            for v = 1:data_stru.data_config.n
                cF_max = max(sim_stru(v).cF);
                Jw_Lp(v) = data_stru.data_config.delP - nRT * theta(3) * cF_max;
            end
            
            if any(Jw_Lp<1e-8)
                trace(p,c) = nan;
                determinant(p,c) = nan;
                min_eig(p,c) = nan;
                cond(p,c) = nan;
                continue
            end
        end    
    end
end

% save results
if ~isfolder('doe_FIM')
    mkdir('doe_FIM');
end
heatmap_doe.cf0 = cf0;
heatmap_doe.delp = delp;
heatmap_doe.trace = trace;
heatmap_doe.determinant = determinant;
heatmap_doe.min_eig = min_eig;
heatmap_doe.cond = cond;
stru_name = ['heatmap_doe-dat', num2str(Dat),'.mat'];
save(stru_name, 'heatmap_doe')
movefile(stru_name,'doe_FIM')

% Visualization
if view_heatmap
    figure
    subplot(1,4,1)
    h1=heatmap(cf0,fliplr(delp),flipud(trace));
    h1.Title = '\fontsize{15} A-optimality';
    h1.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h1.YLabel = '\fontsize{15} Applied pressure [psi]';
    h1.Colormap = parula;
%     h1.ColorbarVisible = 'off';
    h1.MissingDataColor = [0.5, 0.5, 0.5];
    h1.GridVisible = 'off';
    
    subplot(1,4,2)
    h2=heatmap(cf0,fliplr(delp),flipud(determinant));
    h2.Title = '\fontsize{15} D-optimality';
    h2.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h2.YLabel = '\fontsize{15} Applied pressure [psi]';
    h2.Colormap = parula;
%     h2.ColorbarVisible = 'off';
    h2.MissingDataColor = [0.5, 0.5, 0.5];
    h2.GridVisible = 'off';

    subplot(1,4,3)
    h3=heatmap(cf0,fliplr(delp),flipud(min_eig));
    h3.Title = '\fontsize{15} E-optimality';
    h3.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h3.YLabel = '\fontsize{15} Applied pressure [psi]';
    h3.Colormap = parula;
%     h3.ColorbarVisible = 'off';
    h3.MissingDataColor = [0.5, 0.5, 0.5];
    h3.GridVisible = 'off';
    
    subplot(1,4,4)
    h3=heatmap(cf0,fliplr(delp),flipud(cond));
    h3.Title = '\fontsize{15} Modified E-optimality';
    h3.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h3.YLabel = '\fontsize{15} Applied pressure [psi]';
    h3.Colormap = parula;
%     h3.ColorbarVisible = 'off';
    h3.MissingDataColor = [0.5, 0.5, 0.5];
    h3.GridVisible = 'off';

    if length(cf0) <  11
        set(gcf,'Units','Inches','Position',[0.01 0.01 32 5])
    else
        set(gcf,'Units','Inches','Position',[0.01 0.01 64 10])
    end

    plotname = ['doe_heatmap-',model_stru.filenamestr,'.png'];
    saveas(gcf,plotname)
    movefile(plotname,'doe_FIM','f')
end  