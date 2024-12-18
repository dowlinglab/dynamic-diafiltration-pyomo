function weight_stru = weight_config(best_weight, data_stru, model_stru)
% Assemble weight structure based on best_weight options
% 
% Arguments: 
%     best_weight - best weight options
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     weight_stru - weight structure for residuals based on best_weight options
%         default - default weight structure
%         statistical - weight structure from error propagation
%         equal_statistical - equal statistical weight by the number of
%         measurements
%         L1 - weight structure from evaluating L1 norm of Pareto set
%         L2 - weight structure from evaluating L2 norm of Pareto set
%         Linf - weight structure from evaluating Linf norm of Pareto set

% Default weight
raw_weight.m = ones(1,size(data_stru.data_config.M_F0,2));
raw_weight.cp = ones(1,size(data_stru.data_config.C_F0,2));
if isfield(data_stru.data_raw,'cF_exp')
    raw_weight.cr = ones(1,size(data_stru.data_config.C_F0,2));
end

% Assemble weight structure
weight_stru.default = weight_assemble(raw_weight, data_stru);

% Weight from error propagation
if best_weight.statistical
    % 0.01g for mass measurements
    raw_weight.m = ones(1,size(data_stru.data_config.M_F0,2))/0.01;
    % 3% for ICP measurements 
    raw_weight.cp = ones(1,size(data_stru.data_config.C_F0,2))/0.03;
    if isfield(data_stru.data_raw,'cF_exp')
        if ~data_stru.conductivity_cF
            % 3% for ICP measurements 
            raw_weight.cr = ones(1,size(data_stru.data_config.C_F0,2))/0.03;
        else
            % 0.3% for conductivity measurements 
            raw_weight.cr = ones(1,size(data_stru.data_config.C_F0,2))/0.003;
        end
    end
    % Assemble weight structure
    weight_stru.statistical = weight_error_assemble(raw_weight, data_stru);
end

% Equal statistical weight
if best_weight.equal_statistical
    % generate equal weight
    weight_stru.equal_statistical = weight_stru.default;
    
    % loop over vials, adjust equal weight to equal statisitical weight
    for i = 1:data_stru.data_config.n        
        % 0.01g for mass measurements
        weight_stru.equal_statistical(i).m = weight_stru.equal_statistical(i).m/0.01;
        % 3% for ICP measurements
        weight_stru.equal_statistical(i).cp = weight_stru.equal_statistical(i).cp ./(0.03*data_stru.data_raw(i).cV_avg);
        
        if isfield(data_stru.data_raw,'cF_exp')      
            if ~data_stru.conductivity_cF
                % 3% for ICP measurements 
                weight_stru.equal_statistical(i).cr = weight_stru.equal_statistical(i).cr ./(0.03*data_stru.data_raw(i).cF_exp);
            else
                if ~isempty(data_stru.data_raw(i).cF_exp)
                    j = ~isnan(data_stru.data_raw(i).cF_exp(:,1));
                    % 0.3% for conductivity measurements 
                    weight_stru.equal_statistical(i).cr = weight_stru.equal_statistical(i).cr ./(0.003*data_stru.data_raw(i).cF_exp(j));
                end
            end
        end    
    end
end        

% Weight from Pareto set
if best_weight.L1 | best_weight.L2 | best_weight.Linf
    % Add checking for pareto_stru availability
    % if not exist
        % Call pareto function
        % pareto_stru = pareto
    if best_weight.L1
        p = 1;
        % Assemble weight structure from optimal p_norm weight
        weight_stru.w_L1 = calculate_p_norm(p,pareto_stru);
    end

    if best_weight.L2
        p = 2;
        % Assemble weight structure from optimal p_norm weight
        weight_stru.w_L2 = calculate_p_norm(p,pareto_stru);
    end

    if best_weight.Linf
        p = 1e6;
        % Assemble weight structure from optimal p_norm weight
        weight_stru.w_Linf = calculate_p_norm(p,pareto_stru);
    end
end

end

function weight_stru_= weight_p_norm(p,pareto_stru)
    % p-norm distance 
    % Lp_norm = sum(.^p,2).^(1/p);% sum over rows        
    % Find nearest point to Utopia point
    % ind = Lp_norm==min(min(Lp_norm));
    % raw_weight = % optimal weights
    % weight_stru_ = weight_assemble(raw_weight, data_stru);
end