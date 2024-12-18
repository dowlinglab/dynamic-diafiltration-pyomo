function residuals_lsqnonlin = calc_residuals_lsqnonlin(scale_opt, weight_stru_, theta, data_stru, model_stru)
% Calculate weighted residuals needed for fitting with lsqnonlin()
% 
% Arguments: 
%     scale_opt - scaling option 
%     weight_stru_ - assembled weight - substructure in weight_stru
%     theta - fitting parameters
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     residuals_lsqnonlin - weighted residuals needed for fitting with lsqnonlin()

% Calculate unweighted residuals
residuals = calc_residuals(theta, data_stru, model_stru);

r_m = [];
r_cp = [];
r_cr = [];

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

% assemble weighted residuals based on scaling option
switch scale_opt
    case 0        
        
    case 1
        % rescale by Utopia and Alternate Nadir points
        r_m = (1 / (model_stru.fNadirAlt.m - model_stru.fUtopia.m)) .^0.5 .* r_m;
        r_cp = (1 / (model_stru.fNadirAlt.cp - model_stru.fUtopia.cp)) .^0.5 .* r_cp;
        r_cr = (1 / (model_stru.fNadirAlt.cr - model_stru.fUtopia.cr)) .^0.5 .* r_cr ;

end

% assemble residuals to one column
% 100 times to amplify the objectives
residuals_lsqnonlin = 100*[reshape(r_m,[],1);... 
                       reshape(r_cp,[],1);... 
                       reshape(r_cr,[],1)];        

end


