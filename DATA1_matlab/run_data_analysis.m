clear all
%% Load data set and model configuration
% Select dataset
% 501.1 filtration w/ conductivity measurements
% 501.11 filtration w/o conductivity measurements
% 511.12 diafiltration w/ conductivity measurements
% 511.11 diafiltration w/o conductivity measurements

% dat = [501.12];
dat = [501.1,501.11,511.12,511.11];

Lp_w = 3.66; % NF90.5 pure water permeability
mod = '201cvmv'; % default for DATA

% perform default regression
reg_method = 1;
scale_opt = 0;
% Selecting the best weight
best_weight.statistical = false; % Error propagation
best_weight.equal_statistical = true; % Hybrid weight
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm

for i=1:length(dat)
Dat = dat(i)
data_stru = load(['data/data_stru-dataset', num2str(Dat), '.mat']).data_stru;

% create model structure with default options
model_stru = model_config(mod,data_stru);
model_stru.concpolar = true;
model_stru = update_model_config(data_stru,model_stru);

cF = readmatrix([data_stru.filename,'_cF.csv']);
if ~any(Dat==[501.12,511.12])
    model_stru.initialization.C_F0 = cF(1,2:end);
elseif Dat==511.12
    model_stru.initialization.C_F0 = cF(cF(:,1)==380,2:end);
end 
% model_stru = update_model_config_pervial(data_stru,model_stru);


% %% data set preview with initial guess
% % simulation
% simulation.sim_stru = sim_model(model_stru.theta0, data_stru, model_stru);
% plot_sim(data_stru, model_stru, simulation.sim_stru);

%% fitting

weight_stru = weight_config(best_weight, data_stru, model_stru);
fit_stru = fitting(reg_method, scale_opt, weight_stru.equal_statistical, data_stru, model_stru);

plot_sim(data_stru, model_stru, fit_stru.sim_stru);
% Save simulation results
% T = sim_table(data_stru,model_stru,fit_stru);
disp(['Lp = ',num2str(fit_stru.x(1)),' L / m / m / hr / bar',newline,...
    'B = ',num2str(fit_stru.x(2)),' micrometers/s'])
disp(['sigma = ',num2str(fit_stru.x(3)),' dimensionless'])  

disp(['Obj(m) = ',num2str(fit_stru.obj_ind.unscaled.m),newline,...
    'Obj(cv) = ',num2str(fit_stru.obj_ind.unscaled.cp),newline,...
    'Obj(cr) = ',num2str(fit_stru.obj_ind.unscaled.cr)])    
%%
% % contours
% x_var = 'sigma';
% x_ind = 3;
% 
% y_var = 'Lp';
% y_ind = 1;
% grid_density = 50;
% 
% theta = fit_stru.x
% theta(1) = Lp_w;
% contour_stru = calc_contour_2d(data_stru, model_stru, theta, x_var, x_ind, y_var, y_ind, grid_density);
% plot_contour(data_stru, model_stru, contour_stru)
% 
% %%
% x_var = 'B';
% x_ind = 2;
% 
% y_var = 'Lp';
% y_ind = 1;
% grid_density = 50;
% 
% theta = fit_stru.x
% theta(1) = Lp_w;
% contour_stru = calc_contour_2d(data_stru, model_stru, theta, x_var, x_ind, y_var, y_ind, grid_density);
% plot_contour(data_stru, model_stru, contour_stru)

end
%%
% record simulation results
function T=sim_table(data_stru,model_stru,fit_stru)
vial_num = [];
time = [];
cF = [];
mV = [];
cV = [];
cH = [];
for k = 1:length(fit_stru.sim_stru)
    vial_num = [vial_num; k*ones(size(fit_stru.sim_stru(k).time))];
    time = [time; fit_stru.sim_stru(k).time];
    cF = [cF; fit_stru.sim_stru(k).cF];
    mV = [mV; fit_stru.sim_stru(k).mV];
    cV = [cV; fit_stru.sim_stru(k).cV];
    cH = [cH; fit_stru.sim_stru(k).cH];
end

% output simulation result to .csv file
T = table(vial_num,time,cF,mV,cV,cH, 'VariableNames',...
             {'vial_num','time','cF','mV','cV','cH'});
writetable(T,[model_stru.filenamestr,'\fitting-','dat',num2str(data_stru.dataset),'.csv'])  
end
