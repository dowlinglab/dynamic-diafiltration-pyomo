function obj_ind = calc_ind_objectives(scale_opt, theta, data_stru, model_stru)
% Calculate individual objectives
%
% Arguments: 
%     scale_opt - scaling option 
%     weight_stru_ - assembled weight - substructure in weight_stru
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     objective_fmincon - weighted objective needed for fitting with fmincon()

% Calculate unweighted residuals
residuals = calc_residuals(theta, data_stru, model_stru);

best_weight.statistical = false; % Error propagation
best_weight.equal_statistical = false; % Hybrid weight
best_weight.L1 = false; % Minimize L1 norm
best_weight.L2 = false; % Minimize L2 norm
best_weight.Linf = false; % Minimize Linf norm
weight_stru = weight_config(best_weight, data_stru, model_stru);
weight_stru_ = weight_stru.default;

r_m = [];
r_cp = [];
r_cr = [];

% without weight
% % loop over vials, store weighted residuals into array
% for i = 1:data_stru.data_config.n
%     % save mass measurements
%     r_m = [r_m; residuals(i).res_m];
%     % save permeate concentration weighted residuals
%     r_cp = [r_cp; residuals(i).res_cp];
% 
% end
% 
% % if retentate concentration is measured
% if isfield(data_stru.data_raw,'cF_exp')
%     for i = 1:data_stru.data_config.n
%         % if retentate concentration is measured for each vial
%         if ~isempty(data_stru.data_raw(i).cF_exp)
%             % save retentate concentration weighted residuals
%             r_cr = [r_cr;residuals(i).res_cr];
%         end
%     end
% end

% with weight
% loop over vials, store weighted residuals into array
for i = 1:data_stru.data_config.n
    % save mass measurements
    r_m = [r_m; weight_stru_(i).m .* residuals(i).res_m];
    % save permeate concentration weighted residuals
    r_cp = [r_cp; weight_stru_(i).cp .* residuals(i).res_cp];

end

% if retentate concentration is measured
if isfield(data_stru.data_raw,'cF_exp')
    for i = 1:data_stru.data_config.n
        % if retentate concentration is measured for each vial
        if ~isempty(data_stru.data_raw(i).cF_exp)
            % save retentate concentration weighted residuals
            r_cr = [r_cr; weight_stru_(i).cr .* residuals(i).res_cr];
        end
    end
end

obj_ind.unscaled.m = sum(r_m .^2);
obj_ind.unscaled.cp = sum(r_cp .^2);
if isfield(data_stru.data_raw,'cF_exp')
    obj_ind.unscaled.cr = sum(r_cr .^2);
end

% assemble weighted residuals based on scaling option
switch scale_opt
    case 0 
        
    case 1
        % rescale by Utopia and Alternate Nadir points
        obj_ind.scaled.m = (obj_ind.unscaled.m - model_stru.fUtopia.m)./(model_stru.fNadirAlt.m - model_stru.fUtopia.m);
        obj_ind.scaled.cp = (obj_ind.unscaled.cp - model_stru.fUtopia.cp)./(model_stru.fNadirAlt.cp - model_stru.fUtopia.cp);
        obj_ind.scaled.cr = (obj_ind.unscaled.cr - model_stru.fUtopia.cr)./(model_stru.fNadirAlt.cr - model_stru.fUtopia.cr);
end

end

