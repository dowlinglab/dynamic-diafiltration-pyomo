function residuals = calc_residuals(theta, data_stru, model_stru)
% Save simulation residuals into a structure
% 
% Arguments: 
%     theta - fitting parameters
%     data_stru - all experimental data
%     model_stru - partial model configuration
% 
% Returns:
%     residuals - a structure includes different kinds of residuals for
%     each vial
%         res_m - mass residuals
%         res_cp - permeate concentration residuals
%         res_cr - retentate concentration residuals

if ~check_regress_B_per_vial(model_stru)
    % regress overall B
    sim_stru = sim_model(theta, data_stru, model_stru);
else
    % regress B per vial 
    sim_stru = sim_model_pervial(theta, data_stru, model_stru);
end    

% loop over vials, store residuals into a structured array
for i = 1:data_stru.data_config.n
    
    % grab index of reasonable measurements
    j = ~isnan(data_stru.data_raw(i).mass);
    
    % mass residuals
    residuals(i).res_m = data_stru.data_raw(i).mass(j) - sim_stru(i).mV(j);
    % permeate concentration residuals
    residuals(i).res_cp = data_stru.data_raw(i).cV_avg - sim_stru(i).cV(end,:);
    
end

% retentate concentration residuals
% if retentate concentration is measured
if isfield(data_stru.data_raw,'cF_exp')
    for i = 1:data_stru.data_config.n
        % if retentate concentration is measured for each vial
        if ~isempty(data_stru.data_raw(i).cF_exp)
            % one measurement for each vial
            if size(data_stru.data_raw(i).cF_exp,1) == 1
                residuals(i).res_cr = data_stru.data_raw(i).cF_exp - sim_stru(i).cF(end,:);  
            % continous measurements
            else
                % grab index of reasonable measurements
                j = ~isnan(data_stru.data_raw(i).cF_exp(:,1));
                residuals(i).res_cr = data_stru.data_raw(i).cF_exp(j,:) - sim_stru(i).cF(j,:); 
            end
        else
            residuals(i).res_cr = [];
        end
    end
end

end
