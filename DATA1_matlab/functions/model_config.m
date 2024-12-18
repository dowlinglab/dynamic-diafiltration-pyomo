function model_stru = model_config(legacy_model_name,data_stru)
% Save default model configuration into a structure model_stru
% 
% Arguments:
%     legacy_model_name - numerical model name based on options
%     data_stru - all experimental data
% 
% Returns:
%     model_stru - a structrue which includes
%         titlestr - title string for figures
%         filenamestr - file name string for saving files
%         multicomponent - single/multiple component toggle
%         no_crossinter - assumption Bij = Bji = 0 toggle
%         concpolar - consider concentration polarization
%         cp_negligible - permeate concentration negligibility
%         holdup - assumption hold-up volume toggle
%         holdup_part - hold-up volume considered, [under membrane; tube]
%         combine_cvmv - control the integrate variables 
%         Lp_fixed - fix Lp toggle
%         B_fixed - fix B toggle
%         sigma_fixed - fix sigma toggle
%         ch0_fixed - fix initial hold up concentration toggle
%         initialization - model initialization
%             M_F0 - initial feed mass [g]
%             C_F0 - initial feed concentration [mMol/L]
%             M_V0 - initial vial mass [g]
%             C_V0 - initial vial concentration [mMol/L]
%         var_ind - variable indices
%             mF - indices for feed mass
%             cF - indices for feed concentration
%             mV - indices for vial mass
%             cV - indices for vial concentration
%             cH - indices for hold up concentration

model_stru.titlestr = 'Please run update_model_config()';
model_stru.filenamestr = 'default';

if length(data_stru.data_config.C_F0) == 1
    % single/multiple component toggle
    model_stru.multicomponent = false;
else
    model_stru.multicomponent = true;
    % assume Bij = Bji = 0
    model_stru.no_crossinter = true;
end

model_stru.concpolar = false;

% Transfer legacy model name to toggles
switch legacy_model_name
    case{'101cv','101.1cv'}
        % permeate concentration negligibility
        model_stru.cp_negligible = false; 
        % assume a hold-up volume
        model_stru.holdup = false;
        % control the integrate variables 
        model_stru.combine_cvmv = false;
    case{'101cvmv','101.1cvmv'}    
        model_stru.cp_negligible = false; 
        model_stru.holdup = false;
        model_stru.combine_cvmv = true;
    case{'1cv','1.1cv'}
        model_stru.cp_negligible = true; 
        model_stru.holdup = false;
        model_stru.combine_cvmv = false;
    case{'1cvmv','1.1cvmv'}  
        model_stru.cp_negligible = true; 
        model_stru.holdup = false;
        model_stru.combine_cvmv = true;
    case{'2cv','2.1cv'}
        model_stru.cp_negligible = false; 
        model_stru.holdup = true;
        model_stru.holdup_part = [true;true]; %[under membrane; tube]
        model_stru.combine_cvmv = false;
    case{'2cvmv','2.1cvmv'} 
        model_stru.cp_negligible = false; 
        model_stru.holdup = true;
        model_stru.holdup_part = [true;true]; %[under membrane; tube]
        model_stru.combine_cvmv = true;        
    case{'201cv','201.1cv'}
        model_stru.cp_negligible = false; 
        model_stru.holdup = true;
        model_stru.holdup_part = [true;false]; %[under membrane; tube]
        model_stru.combine_cvmv = false;
    case{'201cvmv','201.1cvmv'} 
        model_stru.cp_negligible = false; 
        model_stru.holdup = true;
        model_stru.holdup_part = [true;false]; %[under membrane; tube]
        model_stru.combine_cvmv = true;
    otherwise
        disp('Use default model: 201cvmv')
        model_stru.cp_negligible = false; 
        model_stru.holdup = true;
        model_stru.holdup_part = [true;false]; %[under membrane; tube]
        model_stru.combine_cvmv = true;     
end

% fixed parameter switches
% fix Lp
model_stru.Lp_fixed = false(size(data_stru.data_config.Lp0));  
% fix B
model_stru.B_fixed = ~eye(size(data_stru.data_config.B0));
% fix sigma
model_stru.sigma_fixed = false(size(data_stru.data_config.sigma0)); 
% fix initial hold up concentration
model_stru.ch0_fixed = true(size(data_stru.data_config.C_F0)); 

% initialize model
if data_stru.mode == 'F'
    % filtration mode initialization
    model_stru.initialization.M_F0 = data_stru.data_config.M_F0;
    model_stru.initialization.C_F0 = data_stru.data_config.C_F0;

elseif data_stru.mode == 'D'
    % diafiltration initialization with overflow
    model_stru.initialization.M_F0 = data_stru.data_config.M_F0 + data_stru.data_config.M_O;
    model_stru.initialization.C_F0 = (data_stru.data_config.C_F0 * data_stru.data_config.M_F0...
                                        + data_stru.data_config.C_D * data_stru.data_config.M_O)...
                                        /model_stru.initialization.M_F0;
end
model_stru.initialization.M_V0 = 0.05; % one drop in vial
model_stru.initialization.C_V0 = zeros(size(data_stru.data_raw(1).cV_avg)).';

% update initial retentate concentration in experiments w/ continous conductivity measurements
if size(data_stru.data_raw(1).cF_exp,1) > 1
    model_stru.initialization.C_F0 = data_stru.data_raw(1).cF_exp(1,:);
end

% save indices for variables
% Number of components
nc = data_stru.data_config.nc;
ind(1,:) = [0,0];
% Only filtration mode need mF as an integrated variable, diafiltration
% assume constant mF
if data_stru.mode == 'F'
    ind(1,:) = [1,1];
    model_stru.var_ind.mF = ind(1,:);
end    
ind(2,:)= [ind(1,2) + 1, ind(1,2) + nc];
model_stru.var_ind.cF = ind(2,:);
ind(3,:)= [ind(2,2) + 1, ind(2,2) + 1];
model_stru.var_ind.mV = ind(3,:);
ind(4,:)= [ind(3,2) + 1, ind(3,2) + nc];
model_stru.var_ind.cV = ind(4,:);
% repeat assign indices for cH in update_model_config.m
if model_stru.holdup
    ind(5,:)= [ind(4,2) + 1, ind(4,2) + nc];
    model_stru.var_ind.cH = ind(5,:);
end
end

