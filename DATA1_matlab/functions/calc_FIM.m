function doe_stru = calc_FIM(scale_opt, weight_stru_, theta, step, formula, data_stru, model_stru)
% Calculate FIM using finite difference
% 
% Arguments: 
%     scale_opt - scaling option 
%     weight_stru_ - assembled weight - substructure in weight_stru
%     theta - fitting parameters to evaluate FIM
%     step - step size for finite difference
%     data_stru - all experimental data
%     model_stru - model configuration
% 
% Returns:
%     doe_stru - fitting results which includes
%         x - best fitting parameters
%         resnorm - weighted objective value
%         hessian - Hessian matrix of the optimization problem
%         H_eig - eigen values of the Hessian matrix
%         H_vec - eigen vectors of the Hessian matrix

switch formula
    case 'forward'
        theta_p = [theta, theta * (1+step)];        
    case 'backward'
        theta_p = [theta * (1-step), theta];       
    case 'central'
        theta_p = [theta * (1-step), theta * (1+step)];
end

% Use weighted residuals for lsqnonlin()
residuals_lsqnonlin = @(beta) calc_residuals_lsqnonlin(scale_opt, weight_stru_, beta, data_stru, model_stru);

for i=1:length(theta)
    theta_pi = [theta, theta];
    theta_pi(i,:) = theta_p(i,:);
    Jac(:,i) = (residuals_lsqnonlin(theta_pi(:,2)) - residuals_lsqnonlin(theta_pi(:,1))) ...
            /step;
%     step = 1e-8;
%     theta_1=theta;
%     theta_1(i)=theta_1(i)+step;
%     Jac(:,i) = (residuals_lsqnonlin(theta_1) - residuals_lsqnonlin(theta)) ...
%             /step;
end

if formula == 'central'
    doe_stru.Jac = Jac/2;
else
    doe_stru.Jac = Jac;
end

doe_stru.FIM = Jac' * Jac;
[doe_stru.eig_vec,doe_stru.eig_val] = eig(doe_stru.FIM,'vector');
doe_stru.trace = trace(doe_stru.FIM);
doe_stru.det = det(doe_stru.FIM);
doe_stru.min_eig = min(doe_stru.eig_val);
doe_stru.cond = max(doe_stru.eig_val)/min(doe_stru.eig_val);

end