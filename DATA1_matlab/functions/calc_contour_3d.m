function contour3d_stru = calc_contour_3d(data_stru, model_stru, theta, var, var_ind, grid_density)
% Perform 3D sensitivity analysis-single salt
% 
% Arguments:   
%     data_stru - all experimental data [struc]
%     model_stru - model configuration [struc]
%     theta - point to center sensitivity analysis, the elements must
%     be the same as model_stru.theta
%     var - name of variables to perturb [str], this must be a
%     field in theta
%     var_ind - variable position indices in theta
%
% Returns:
%     contour3d_stru - structure save all contours data
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

if ~model_stru.multicomponent
    var = {'B';'Lp';'sigma'};
    var_ind = [2;1;3];
end

contour3d_stru.var = var;
contour3d_stru.grid_density = grid_density;

% Generate nd grid
var_lb = model_stru.theta_lb(var_ind);
var_ub = model_stru.theta_ub(var_ind);

% update range for Lp
var_lb(var_ind(var == "Lp")) = 0.1 * theta(var_ind(var == "Lp"));
var_ub(var_ind(var == "Lp")) = 2 * theta(var_ind(var == "Lp"));

% update range for B
var_lb(var_ind(var == "B")) = 0;
var_ub(var_ind(var == "B")) = 2;

[X,Y,Z] = ndgrid(linspace(var_lb(var_ind(var == "B")),var_ub(var_ind(var == "B")),grid_density),...
    linspace(var_lb(var_ind(var == "Lp")),var_ub(var_ind(var == "Lp")),grid_density),...
    linspace(var_lb(var_ind(var == "sigma")),var_ub(var_ind(var == "sigma")),6));

% Allocate arrays for contour plot
xx = X(:); yy = Y(:); zz = Z(:);
f_m = zeros(length(xx),data_stru.data_config.nc);
f_cp = zeros(length(xx),data_stru.data_config.nc);
f_cr = zeros(length(xx),data_stru.data_config.nc);
contour3d_stru.xx = xx;
contour3d_stru.yy = yy;
contour3d_stru.zz = zz;

beta = theta;
% Loop over parameter value combinations
for i=1:length(xx)
    beta(var_ind(1)) = xx(i);
    beta(var_ind(2)) = yy(i);
    beta(var_ind(3)) = zz(i);
    
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
    contour3d_stru.Objectives(i).mass = f_m;
    contour3d_stru.Objectives(i).m_min_ind = f_m==min(min(f_m));
    
    contour3d_stru.Objectives(i).perm_conc = f_cp(:,i);
    contour3d_stru.Objectives(i).cp_min_ind = f_cp(:,i)==min(min(f_cp(:,i)));
    
    if isfield(data_stru.data_raw,'cF_exp')
        contour3d_stru.Objectives(i).reten_conc = f_cr(:,i);
        contour3d_stru.Objectives(i).cr_min_ind = f_cr(:,i)==min(min(f_cr(:,i)));
    end
end

save([model_stru.filenamestr,'\contour3d_stru','.mat'], 'contour3d_stru')

% save objective values in csv file
contourdata = table(xx,yy,zz,f_m,f_cp,f_cr,'VariableNames',...
    {var{1},var{2},var{3},'Obj_mass','Obj_perm_conc','Obj_reten_conc'});
writetable(contourdata,[model_stru.filenamestr,'\contour3ddata-',model_stru.filenamestr,'.csv'])
end

