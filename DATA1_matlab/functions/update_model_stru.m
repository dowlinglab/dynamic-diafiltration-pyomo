function model_stru = update_model_stru(theta, model_stru)
% Use theta to update the value of the parameter
% 
% Arguments:
%     model_stru - model configuration based on default options
% 
% Returns:
%     model_stru - model structrue update with
%         model_para - model parameters
%             Lp - Lp for simulation [L/m^2/bar/h]
%             B - B for simulation [um/s]
%             sigma - sigma for simulation [-]
%             chi - initial hold up concentration for simulation [mMol/L]


ind = model_stru.para_ind;

% unpack fitted parameters
% select parameters to fit from .._fixed toggles
model_stru.model_para.Lp(~model_stru.Lp_fixed) = theta(ind.Lp(1):ind.Lp(2));

model_stru.model_para.B(~model_stru.B_fixed) = theta(ind.B(1):ind.B(2));
% reshape model_para.B in matrix form
model_stru.model_para.B = reshape(model_stru.model_para.B,size(model_stru.B0));

model_stru.model_para.sigma(~model_stru.sigma_fixed) = theta(ind.sigma(1):ind.sigma(2));

if ~all(model_stru.ch0_fixed)
    model_stru.model_para.chi(~model_stru.ch0_fixed) = theta(ind.chi(1):ind.chi(2));
end

end

