function contour_stru = calc_contour_2d(data_stru, model_stru, theta, x_var, x_ind, y_var, y_ind, grid_density)
% Perform 2D sensitivity analysis
% 
% Arguments:   
%     data_stru - all experimental data [struc]
%     model_stru - model configuration [struc]
%     theta - point to center sensitivity analysis, the elements must
%     be the same as model_stru.theta
%     x_var - name of first variable to perturb [str], this must be a
%     field in theta
%     x_ind - x_var position index in theta
% %     x_values - values for the first variable [vector]
%     y_var - name of second variable to perturb [str], this must be a
%     field in theta
%     y_ind - y_var position index in theta
% %     y_values - values for the second variable [vector]
%
% Returns:
%     contour_stru - structure save all contours data
%         x_var - name of first variable for sensitivty analysis
%         y_var - name of second variable for sensitivity analysis
%         xx - values for x_var in sensitivity analysis [vector]
%         yy - values for y_var in sensitvity analysis [vector]
%         Objectives - contains arrary structure
%             mass - log10 transform least square mass residual objective values[vector]
%             m_min_ind - minimum of mass objective
%             perm_conc - log10 transform least square permeate concentration residual objective[vector]
%             cp_min_ind - minimum of permeate conccentration objective
%             reten_conc - log10 transform least square retentate concentration residual objective[vector]
%             cr_min_ind - minimum of retentate conccentration objective

contour_stru.x_var = x_var;
contour_stru.y_var = y_var;
contour_stru.grid_density = grid_density;

% Generate mesh grid
x_lb = model_stru.theta_lb(x_ind);
y_lb = model_stru.theta_lb(y_ind);
x_ub = model_stru.theta_ub(x_ind);
y_ub = model_stru.theta_ub(y_ind);
% update range for Lp
if x_var == "Lp"
    x_lb = 0.1 * theta(x_ind);
    x_ub = 2 * theta(x_ind);
end

if y_var == "Lp"
    y_lb = 0.1 * theta(y_ind);
    y_ub = 2 * theta(y_ind);
end

[X,Y] = meshgrid(linspace(x_lb,x_ub,grid_density),...
    linspace(y_lb,y_ub,grid_density));

% Allocate arrays for contour plot
xx = X(:); yy = Y(:);
f_m = zeros(length(xx),data_stru.data_config.nc);
f_cp = zeros(length(xx),data_stru.data_config.nc);
f_cr = zeros(length(xx),data_stru.data_config.nc);
contour_stru.xx = xx;
contour_stru.yy = yy;

beta = theta;
% Loop over parameter value combinations
for i=1:length(xx)
    beta(x_ind) = xx(i);
    beta(y_ind) = yy(i);

    % Evaluate objectives
    scale_opt = 0; %no scaling
    obj_ind = calc_ind_objectives(scale_opt, beta, data_stru, model_stru);
    
    % save log10 transformed individual objectives
    f_m(i,:) = log10(obj_ind.unscaled.m);
    f_cp(i,:) = log10(obj_ind.unscaled.cp);
    if isfield(data_stru.data_raw,'cF_exp')
        f_cr(i,:) = log10(obj_ind.unscaled.cr);  
    end
end

% save objective values
% find index corresponds minimum on the contour plot
% loop over components
for i=1:data_stru.data_config.nc
    contour_stru.Objectives(i).mass = f_m;
    contour_stru.Objectives(i).m_min_ind = f_m==min(min(f_m));
    
    contour_stru.Objectives(i).perm_conc = f_cp(:,i);
    contour_stru.Objectives(i).cp_min_ind = f_cp(:,i)==min(min(f_cp(:,i)));
    
    if isfield(data_stru.data_raw,'cF_exp')
        contour_stru.Objectives(i).reten_conc = f_cr(:,i);
        contour_stru.Objectives(i).cr_min_ind = f_cr(:,i)==min(min(f_cr(:,i)));
    end
end

% save([model_stru.filenamestr,'\contour_stru-','x_',x_var,'-y_',y_var, '.mat'], 'contour_stru')

% save objective values in csv file
contourdata = table(xx,yy,f_m,f_cp,f_cr, 'VariableNames',...
    {x_var,y_var,'Obj_mass','Obj_concentration','Obj_retentate_concentration'});
writetable(contourdata,[model_stru.filenamestr,'\contourdata-','x_',x_var,'-y_',y_var,'.csv'])
end

