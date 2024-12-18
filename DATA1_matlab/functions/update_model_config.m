function model_stru = update_model_config(data_stru,model_stru)
% Update model configuration based on customized model options
%
% Arguments:
%     data_stru - all experimental data
%     model_stru - model configuration based on default options
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
%             C_H0 - initial hold up concentration [mMol/L]
%         var_ind - variable indices
%             mF - indices for feed mass
%             cF - indices for feed concentration
%             mV - indices for vial mass
%             cV - indices for vial concentration
%             cH - indices for hold up concentration
%         para_ind - parameter indices
%             Lp - indices for Lp value in theta
%             B - indices for B value in theta
%             sigma - indices for sigma value in theta
%             chi - indices for chi value in theta
%         Lp0 - updated initial guess for Lp [L/m^2/bar/h]
%         B0 - updated initial guess for B [um/s]
%         sigma0 - updated initial guess for sigma [-]
%         theta0 - updated target parameters initial guess
%         theta_lb - target parameters lower bound
%         theta_ub - target parameters upper bound
%         model_para - model parameters
%             Lp - Lp for simulation [L/m^2/bar/h]
%             B - B for simulation [um/s]
%             sigma - sigma for simulation [-]
%             chi - initial hold up concentration for simulation [mMol/L]

% Specify dataset, model and integrate variables for titles and file names
model_stru.titlestr = ['dataset ', num2str(data_stru.dataset)];
model_stru.filenamestr = ['dat ', num2str(data_stru.dataset)];

if length(data_stru.data_config.C_F0) == 1
    % single/multiple component toggle
    model_stru.multicomponent = false;
    model_stru.model_name = 'oneCPNT';
    model_stru.titlestr = [model_stru.titlestr, ' single component'];
    model_stru.filenamestr = [model_stru.filenamestr, ' oneCPNT'];
else 
    model_stru.multicomponent = true;
    model_stru.model_name = 'multiCPNT';
    model_stru.titlestr = [model_stru.titlestr, ' multi-component'];
    model_stru.filenamestr = [model_stru.filenamestr, ' multiCPNT'];
end

if data_stru.mode == 'F'
    model_stru.titlestr = [model_stru.titlestr, ' filtration']; 
elseif data_stru.mode == 'D'
    model_stru.titlestr = [model_stru.titlestr, ' diafiltration'];
end
    
if model_stru.cp_negligible   
    model_stru.model_name = [model_stru.model_name, ' cp_neg'];
    model_stru.titlestr = [model_stru.titlestr, ' assuming c_p negligible'];
    model_stru.filenamestr = [model_stru.filenamestr, ' cp_neg'];   
end    

if model_stru.holdup
    % Holdup mass 
    % Ouimet, J: email on 20-Dec-2019, 1:40pm. Holdup volume = 0.25 mL (under membrane) +
    % 0.8 mL (tube)
    mH = [0.25, 0.8];%[under membrane, tube]
    model_stru.initialization.M_H = mH * model_stru.holdup_part; %[g]
    model_stru.initialization.C_H0 = data_stru.data_raw(1).cV_avg.';
    % save indices for cH as integrated variable
    model_stru.var_ind.cH = [model_stru.var_ind.cV(2) + 1, model_stru.var_ind.cV(2) + data_stru.data_config.nc];
    model_stru.model_name = [model_stru.model_name, ' holdup'];
    model_stru.titlestr = [model_stru.titlestr, ' w/ holdup'];
    model_stru.filenamestr = [model_stru.filenamestr, ' holdup'];
    
    if data_stru.mode == 'D'
        % update diafiltration initialization with overflow in hold up volume
        model_stru.initialization.M_F0 = model_stru.initialization.M_F0;
        model_stru.initialization.C_F0 = (model_stru.initialization.C_F0 * model_stru.initialization.M_F0...
                                            + (data_stru.data_config.C_D - model_stru.initialization.C_H0) * sum(mH))...
                                            /model_stru.initialization.M_F0;
    end
else
    model_stru.initialization.C_H0 = NaN;
end

% interrupt if the model selected is unimplemented
switch model_stru.model_name
    case{'oneCPNT'}
    case{'oneCPNT cp_neg'}
    case{'oneCPNT holdup'}
    case{'multiCPNT'}
    case{'multiCPNT cp_neg'}
    case{'multiCPNT holdup'}
    otherwise
        error('Selected model unimplemented')
end

if model_stru.multicomponent
    if model_stru.no_crossinter
        model_stru.model_name = [model_stru.model_name, ' no_crossinter'];
        model_stru.titlestr = [model_stru.titlestr, ' assuming no cross interaction'];
        model_stru.filenamestr = [model_stru.filenamestr, ' no_crossinter'];
    end
end

if model_stru.concpolar
    model_stru.model_name = [model_stru.model_name, ' concpolar'];
    model_stru.titlestr = [model_stru.titlestr, ' considering concentration polarization'];
    model_stru.filenamestr = [model_stru.filenamestr, ' concpolar'];
    % save indices for cIn as a variable
    model_stru.var_ind.cH = [model_stru.var_ind.cV(2) + 1, model_stru.var_ind.cV(2) + data_stru.data_config.nc];
    if ~model_stru.holdup
        model_stru.var_ind.cIn = [model_stru.var_ind.cV(2) + 1, model_stru.var_ind.cV(2) + data_stru.data_config.nc];
    else
        model_stru.var_ind.cIn = [model_stru.var_ind.cH(2) + 1, model_stru.var_ind.cH(2) + data_stru.data_config.nc];  
    end
end

% append titles and file names with other options
if model_stru.combine_cvmv
    model_stru.model_name = [model_stru.model_name,' cvmv'];
    model_stru.titlestr = [model_stru.titlestr, ' combined variable c_vm_v'];
    model_stru.filenamestr = [model_stru.filenamestr, ' cvmv'];
else
    model_stru.model_name = [model_stru.model_name,' cv'];
    model_stru.titlestr = [model_stru.titlestr, ' single variable c_v'];
    model_stru.filenamestr = [model_stru.filenamestr, ' cv'];
end

if any(model_stru.Lp_fixed) 
    model_stru.model_name = [model_stru.model_name,' fixed Lp'];
    model_stru.titlestr = [model_stru.titlestr, ' fixed Lp'];
    model_stru.filenamestr = [model_stru.filenamestr, ' fixed Lp'];
end
    data_stru.data_config.theta0 = data_stru.data_config.Lp0(~model_stru.Lp_fixed);

if any(diag(model_stru.B_fixed))
    model_stru.model_name = [model_stru.model_name,' fixed B'];
    model_stru.titlestr = [model_stru.titlestr, ' fixed B'];
    model_stru.filenamestr = [model_stru.filenamestr, ' fixed B'];
end

if any(model_stru.sigma_fixed) 
    model_stru.model_name = [model_stru.model_name,' fixed sigma'];
    model_stru.titlestr = [model_stru.titlestr, ' fixed sigma'];
    model_stru.filenamestr = [model_stru.filenamestr, ' fixed sigma'];
end

if ~all(model_stru.ch0_fixed)
    if ~model_stru.holdup
        error('Improper fitting option')
    end
else
    if model_stru.holdup
        model_stru.model_name = [model_stru.model_name,' fixed ch0'];
        model_stru.titlestr = [model_stru.titlestr, ' fixed c_h0'];
        model_stru.filenamestr = [model_stru.filenamestr, ' fixed ch0'];
    end
end

% Parameters initial guess from data library(updated with classical calculation)
% Only take in non-fixed parameters
data_stru.data_config.theta0 = [data_stru.data_config.Lp0(~model_stru.Lp_fixed);...
                                data_stru.data_config.B0(~model_stru.B_fixed);...
                                data_stru.data_config.sigma0(~model_stru.sigma_fixed);...
                                model_stru.initialization.C_H0(~model_stru.ch0_fixed)];

                            
% Create a folder to store files if it doesn't exist
if ~isfolder(model_stru.filenamestr)
    mkdir(model_stru.filenamestr);
end

% if there is a customized initial guess file, update initial guess
if exist([model_stru.filenamestr,'/initial_guess.csv'], 'file')
    T = readtable([model_stru.filenamestr,'/initial_guess.csv']);
    % Update initial guess of parameters
    % One Lp value for multicomponent
    model_stru.Lp0 = T.initialLp(~isnan(T.initialLp)); 
    % Beginning index of B in table
    Bidx = find(string(T.Properties.VariableNames) == "initialB");
    % read in B value for single component/B matrix for multicomponent
    model_stru.B0 = table2array(T(:,Bidx:(Bidx + data_stru.data_config.nc -1)));
    model_stru.sigma0 = T.initialSigma; 
else
    model_stru.Lp0 = data_stru.data_config.Lp0; 
    model_stru.B0 = data_stru.data_config.B0;    
    model_stru.sigma0 = data_stru.data_config.sigma0; 
end

% parameter order: Lp, B(in column), sigma, ch0
% select parameters to fit from .._fixed toggles
B0 = model_stru.B0(~model_stru.B_fixed);
B0 = reshape(B0,[],1);

model_stru.theta0 = [model_stru.Lp0(~model_stru.Lp_fixed);...
                     B0;...
                     model_stru.sigma0(~model_stru.sigma_fixed);...
                     model_stru.initialization.C_H0(~model_stru.ch0_fixed)]; 
% theta lower bound and upper bound
model_stru.theta_lb = [0.1*ones(size(model_stru.Lp0(~model_stru.Lp_fixed)));...
                       1e-6*ones(size(B0));...
                       zeros(size(model_stru.sigma0(~model_stru.sigma_fixed)));...
                       zeros(size(model_stru.initialization.C_H0(~model_stru.ch0_fixed)))];
model_stru.theta_ub = [10*ones(size(model_stru.Lp0(~model_stru.Lp_fixed)));...
                       2*ones(size(B0));...
                       ones(size(model_stru.sigma0(~model_stru.sigma_fixed)));...
                       model_stru.initialization.C_H0(~model_stru.ch0_fixed)];

% initialize model parameters
model_stru.model_para.Lp = model_stru.Lp0; 
model_stru.model_para.B = model_stru.B0;    
model_stru.model_para.sigma = model_stru.sigma0; 

% update initial holdup concentration
if model_stru.holdup   
%     % Applied pressure
%     delP = data_stru.data_config.delP; %[bar]
%     % Gas constant
%     R = 8.314e-5; %[cm^3 bar / micromol / K]
%     % Temperature
%     T = data_stru.data_config.Temp; %[K]
%     % Membrane area
%     Am = data_stru.data_config.Am; %[cm^2]
%     % Density
%     rho = data_stru.data_config.rho; %[g/cm^3]
%     % Number of dissolved species
%     ni = data_stru.data_config.ni;
%     % (L / m / m / hr / bar) * (1 hr / 3600 s) * (1 m / 100 cm) * (1000 cm^3 / L)  
%     Lp = model_stru.model_para.Lp/36000; % cm / bar / s
%     % (micrometer / s) * (1 cm / 1E4 micrometer)
%     B = model_stru.model_para.B/10000; % cm / s
%     sigma = model_stru.model_para.sigma; % dimensionless 
%     cF = model_stru.initialization.C_F0;
%     model_stru.initialization.C_H0 = B*cF /(Lp*(delP - sum(cF.*(ni.*sigma)*R*T)) + B);
    model_stru.initialization.C_H0 = model_stru.initialization.C_H0*0.8;
end

if ~all(model_stru.ch0_fixed)

    model_stru.model_para.chi = model_stru.initialization.C_H0;
end

% save indices for parameters in theta from .._fixed toggles
% Lp
% Assume the size of Lp is unknown but from user set
ind.Lp(1) = 1;
ind.Lp(2) = ind.Lp(1) + sum(~model_stru.Lp_fixed(:) == 1) - 1;
model_stru.para_ind.Lp = ind.Lp;
% B
ind.B(1) = ind.Lp(2) + 1;
ind.B(2) = ind.B(1) + sum(~model_stru.B_fixed(:) == 1) - 1;
model_stru.para_ind.B = ind.B;
% sigma
ind.sigma(1) = ind.B(2) + 1;
ind.sigma(2) = ind.sigma(1) + sum(~model_stru.sigma_fixed(:) == 1) - 1;
model_stru.para_ind.sigma = ind.sigma;
% initial hold up concentration chi
if ~all(model_stru.ch0_fixed)
    ind.chi(1) = ind.sigma(2) + 1;
    ind.chi(2) = ind.chi(1) + sum(~model_stru.ch0_fixed(:) == 1) - 1;
    model_stru.para_ind.chi = ind.chi;
end

% save(['model_stru-Dataset', num2str(Dat), 'model', num2str(model), '.mat', 'data_stru')
end

