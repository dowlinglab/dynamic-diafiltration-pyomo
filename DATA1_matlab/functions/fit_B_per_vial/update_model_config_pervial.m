function model_stru = update_model_config_pervial(data_stru,model_stru)
% Update model configuration for regressing B per vial
%     para_pervial_ind - parameter indices (B for all vials included)
% Rearrange 
%     theta0 - target parameters initial guess 
%     theta_lb - target parameters lower bound
%     theta_ub - target parameters upper bound
% 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - model configuration

% Update file name and folder name
oldfolder = model_stru.filenamestr;
model_stru.filenamestr = [model_stru.filenamestr, ' B_pervial']; 
movefile(oldfolder,model_stru.filenamestr)   

ind = model_stru.para_ind;
% Lp remains the same

% Update indeces for B per vial
% extend B size by number of vials
% size of B per vial
sizeB = ind.B(2) - ind.B(1) + 1;

% Number of vials
n = data_stru.data_config.n; 
% Loop over vials
for i = 2:n 
    % Store begin and end indices of B for ith vial in ith row of para_ind.B    
    ind.B(i,1) = ind.B(i-1,1) + sizeB;
    ind.B(i,2) = ind.B(i-1,2) + sizeB;
end

% sigma
sizeSigma = ind.sigma(2) - ind.sigma(1) + 1;
ind.sigma(1) = ind.B(end,2) + 1;
ind.sigma(2) = ind.sigma(1) + sizeSigma - 1;
model_stru.para_pervial_ind.sigma = ind.sigma;

% initial hold up concentration chi
if ~all(model_stru.ch0_fixed)
    sizechi = ind.chi(2) - ind.chi(1) + 1; 
    ind.chi(1) = ind.sigma(2) + 1;
    ind.chi(2) = ind.chi(1) + sizechi - 1;
end

% Add parameter indices into model_stru
model_stru.para_pervial_ind = ind;

% Rearrange theta0, theta_lb, theta_ub
% parameter order: Lp, B(in column for each vial), sigma, ch0
ind = model_stru.para_ind;

theta0 = model_stru.theta0;
theta_lb = model_stru.theta_lb;
theta_ub = model_stru.theta_ub;

model_stru.theta0 = [theta0(ind.Lp(1):ind.Lp(2));...
        repmat(theta0(ind.B(1):ind.B(2)),n,1);... % Repeat n copies of B in theta0
        theta0(ind.sigma(1):end)];

model_stru.theta_lb = [theta_lb(ind.Lp(1):ind.Lp(2));...
        repmat(theta_lb(ind.B(1):ind.B(2)),n,1);... % Repeat n copies of B in theta_lb
        theta_lb(ind.sigma(1):end)];    
    
model_stru.theta_ub = [theta_ub(ind.Lp(1):ind.Lp(2));...
        repmat(theta_ub(ind.B(1):ind.B(2)),n,1);... % Repeat n copies of B in theta_ub
        theta_ub(ind.sigma(1):end)];    

end