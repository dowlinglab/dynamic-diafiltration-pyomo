function fit_stru = fitting(reg_method, scale_opt, weight_stru_, data_stru, model_stru)
% Perform regression
% 
% Arguments: 
%     reg_method - regression method [integer]
%                   1:lsqnonlin()
%                   2:fmincon()
%     scale_opt - scaling option 
%     weight_stru_ - assembled weight - substructure in weight_stru
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     fit_stru - fitting results which includes
%         x - best fitting parameters
%         resnorm - weighted objective value
%         hessian - Hessian matrix of the optimization problem
%         H_eig - eigen values of the Hessian matrix
%         H_vec - eigen vectors of the Hessian matrix

% lower and upper bounds for parameters
% parameter order: Lp, B, sigma
lb = model_stru.theta_lb;
ub = model_stru.theta_ub;

% Perform regression
if reg_method == 1
    % Use lsqnonlin()
%     option.FiniteDifferenceStepSize = 1e-8
    residuals_lsqnonlin = @(beta) calc_residuals_lsqnonlin(scale_opt, weight_stru_, beta, data_stru, model_stru);
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqnonlin(residuals_lsqnonlin,model_stru.theta0,lb,ub);
    fit_stru.objective = resnorm;
    fit_stru.jacobian = full(jacobian);
    fit_stru.FIM = full(jacobian' * jacobian);
    
else reg_method == 2
    % Use fmincon()
    objective_fmincon = @(beta) calc_objective_fmincon(scale_opt, weight_stru_, beta, data_stru, model_stru);
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective_fmincon,model_stru.theta0,[],[],[],[],lb,ub,[]);
    %hessian
    %[H_vec,H_eig] = eig(hessian,'vector')
end
        
fit_stru.x = x;

% unpack fitted parameters
% select parameters to fit from .._fixed toggles
if ~check_regress_B_per_vial(model_stru)
    % regress overall B
    ind = model_stru.para_ind;
    B = model_stru.B0;
    B(~model_stru.B_fixed) = x(ind.B(1):ind.B(2));
    % reshape model_para.B in matrix form
    fit_stru.B = reshape(B,size(model_stru.B0));
else
    % regress B per vial
    ind = model_stru.para_pervial_ind;    
    for i = 1:data_stru.data_config.n
        B = model_stru.B0;
        B(~model_stru.B_fixed) = x(ind.B(i,1):ind.B(i,2));
        % reshape model_para.B in matrix form
        fit_stru.B{i} = reshape(B,size(model_stru.B0));     
    end    
end    
    
fit_stru.Lp(~model_stru.Lp_fixed) = x(ind.Lp(1):ind.Lp(2));

fit_stru.sigma(~model_stru.sigma_fixed) = x(ind.sigma(1):ind.sigma(2));

if ~all(model_stru.ch0_fixed)
    fit_stru.chi(~model_stru.ch0_fixed) = x(ind.chi(1):ind.chi(2));
end

% Save simulation resulta
if ~check_regress_B_per_vial(model_stru)
    % regress overall B
    fit_stru.sim_stru = sim_model(x, data_stru, model_stru);
else
    % regress B per vial 
    fit_stru.sim_stru = sim_model_pervial(x, data_stru, model_stru);
end

fit_stru.obj_ind = calc_ind_objectives(scale_opt, x, data_stru, model_stru);

save([model_stru.filenamestr,'\fit_stru','.mat'], 'fit_stru')
end




