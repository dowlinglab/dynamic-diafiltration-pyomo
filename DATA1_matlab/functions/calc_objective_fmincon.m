function objective_fmincon = calc_objective_fmincon(scale_opt, weight_stru_, theta, data_stru, model_stru)
% Calculate objective needed for fitting with fmincon()
% 
% Arguments: 
%     scale_opt - scaling option 
%     weight_stru_ - assembled weight - substructure in weight_stru
%     theta - fitting parameters
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     objective_fmincon - weighted objective needed for fitting with fmincon()

% Calculate weighted residuals
residuals_lsqnonlin = calc_residuals_lsqnonlin(scale_opt, weight_stru_, theta, data_stru, model_stru);
% assemble least square objective
objective_fmincon = sum(residuals_lsqnonlin .^2);

end

