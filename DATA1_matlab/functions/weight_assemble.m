function weight_stru_ = weight_assemble(raw_weight, data_stru)
% Assemble default weight structure for residuals
% 
% Arguments: 
%     raw_weight - scalar weights of each objective
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     weight_stru_ - weights correspond each residual

% number of cr measurements
n_crmeas = 0;

% loop over vials, store weights into a structured array
for i = 1:data_stru.data_config.n
    
    % grab index of reasonable measurements
    j = ~isnan(data_stru.data_raw(i).mass);
    
    % weights of mass residuals (average weight - per vial per data)
    weight_stru_(i).m = sqrt(raw_weight.m * ones(size(data_stru.data_raw(i).mass(j)))/size(data_stru.data_raw(i).mass(j),1)/data_stru.data_config.n);
    % weights of permeate concentration residuals (average weight - 1 data per via1)
    weight_stru_(i).cp = sqrt(raw_weight.cp .* ones(size(data_stru.data_raw(i).cV_avg))/data_stru.data_config.n);

    % if retentate concentration is measured
    if isfield(data_stru.data_raw,'cF_exp')
        if ~isempty(data_stru.data_raw(i).cF_exp)
            % count vial with cr measurements 
            n_crmeas = n_crmeas + 1; 
        end
    end
end

% weights of retentate concentration residuals (average weight - counted cr measurements)
% if retentate concentration is measured
if isfield(data_stru.data_raw,'cF_exp')
    for i = 1:data_stru.data_config.n
        % if retentate concentration is measured for each vial
        if ~isempty(data_stru.data_raw(i).cF_exp)
            % one measurement for each vial
            if size(data_stru.data_raw(i).cF_exp,1) == 1
                weight_stru_(i).cr = sqrt(raw_weight.cr .* ones(size(data_stru.data_raw(i).cV_avg))/n_crmeas);
            % continous measurements
            else
                % grab index of reasonable measurements
                j = ~isnan(data_stru.data_raw(i).cF_exp(:,1));
                weight_stru_(i).cr = sqrt(raw_weight.cr .* ones(size(data_stru.data_raw(i).cF_exp(j)))/size(data_stru.data_raw(i).cF_exp(j),1)/n_crmeas);
            end
        else
            weight_stru_(i).cr = [];
        end
    end
end

end

