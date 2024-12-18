function rhs = diafiltration_ode(t, z, data_stru, model_stru)
% Include odes for all model
% 
% Arguments:
%     t - simulation time
%     z - integrate variables
%     data_stru - all experimental data
%     model_stru - model configuration based on default options

% 
% Returns:
%     rhs - left hand side of odes

% units for parameter estimation
% B: micrometers / s
% Lp: L / m / m / hr / bar
% sigma: dimensionless

% units for ODE model
% B: cm /  s
% Lp: cm / bar / s
% sigma: dimensionless

% Applied pressure
delP = data_stru.data_config.delP; %[bar]

% Gas constant
R = 8.314e-5; %[cm^3 bar / micromol / K]

% Temperature
T = data_stru.data_config.Temp; %[K]

% Membrane area
Am = data_stru.data_config.Am; %[cm^2]

% Density
rho = data_stru.data_config.rho; %[g/cm^3]

% Kinematic viscosity
nu = 8.927e-3; %[cm^2/s]

% Diameter of stirred cell
% b = 2.54; %[cm]
b = 2.2860;

% Average velocity within the system
% Add to data library
v = 350/60 * pi * b ; %[cm/s]

% Diffusion coefficient
%[cm^2/s]
for i = 1:size(data_stru.data_config.namec,1)
    if data_stru.data_config.namec(i,1)=='K'
        D(i,:) = 1.960e-5;
    elseif data_stru.data_config.namec(i,1)=='M'
        D(i,:) = 0.705e-5;
    end
end

% Number of dissolved species
ni = data_stru.data_config.ni;

% unpack fitted parameters and convert to units for ODE model
% (L / m / m / hr / bar) * (1 hr / 3600 s) * (1 m / 100 cm) * (1000 cm^3 / L)  
Lp = model_stru.model_para.Lp/36000; % cm / bar / s
% (micrometer / s) * (1 cm / 1E4 micrometer)
B = model_stru.model_para.B/10000; % cm / s
sigma = model_stru.model_para.sigma; % dimensionless 

% unpack variables
ind = model_stru.var_ind;
% filtration mode
if data_stru.mode == 'F'
    mF = z(ind.mF(1):ind.mF(2));
% diafiltration mode
elseif data_stru.mode == 'D'
    % constant mF
    mF = model_stru.initialization.M_F0;
end 

cF = z(ind.cF(1):ind.cF(2));

mV = z(ind.mV(1):ind.mV(2));

if model_stru.combine_cvmv
    % may cause the issue when switch vials
    % (mv*cv) /mv
    cV = z(ind.cV(1):ind.cV(2))./mV;
else
    cV = z(ind.cV(1):ind.cV(2));
end

% create RHS
rhs = zeros(size(z));

% model without holdup
if ~model_stru.holdup    
    % intermediate
    % [cm / bar / s] * ([bar] - [cm^3 bar / micromol / K] * [K] * [micromol / cm^3])
    % = [cm / s]
    if model_stru.cp_negligible  
        Jw = Lp*(delP - sum(cF.*(ni.*sigma)*R*T));
    else
        Jw = Lp*(delP - sum((cF - cV).*(ni.*sigma)*R*T));
    end

    % intermediate
    % [cm / s] * [micromol / cm^2]
    % = [micromol / cm^2 / s]
    Js = B*(cF - cV);  
    
    % consider concentration polarization
    if model_stru.concpolar
        cIn = z(ind.cIn(1):ind.cIn(2));
        
        % intermediate
        % [cm / bar / s] * ([bar] - [cm^3 bar / micromol / K] * [K] * [micromol / cm^3])
        % = [cm / s]
        Jw = Lp*(delP - sum((cIn - cV).*(ni.*sigma)*R*T));

        % intermediate
        % [cm / s] * [micromol / cm^2]
        % = [micromol / cm^2 / s]
        Js = B*(cIn - cV);   
        
        k = 0.23 * v^0.57 .* D.^0.67 / (nu^0.24 * b^0.43);               
        rhs(ind.cIn(1):ind.cIn(2)) = (cF - cV) .* exp(Jw./ k) + cV - cIn;   
    end

    % dcF_dt
    % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^3] * [cm / s] - [micromol / cm^2 / s])
    % = [micromol / cm^3 / s] 
    if data_stru.mode == 'F'
        rhs(ind.cF(1):ind.cF(2)) = Am * rho / mF * (cF * Jw - Js);
    elseif data_stru.mode == 'D'
        % Dialysate concentration
        cD = data_stru.data_config.C_D;
        rhs(ind.cF(1):ind.cF(2)) = Am * rho / mF * (cD * Jw - Js);
    end    

    if model_stru.combine_cvmv
        % dmVcV_dt = JsAmrho
        % [micromol / cm^2 / s] * [cm^2] * [g / cm^3]
        % = [micromol * g / cm^3 / s]      
        rhs(ind.cV(1):ind.cV(2)) = Js * Am * rho;
    else
        % dcV_dt
        % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^2 / s] - [micromol / cm^3] * [cm / s]) 
        % = [micromol / cm^3 / s]
        rhs(ind.cV(1):ind.cV(2)) = Am * rho / mV * (Js - Jw*cV);
    end

% model with holdup
elseif model_stru.holdup

    mH = model_stru.initialization.M_H; %[g] 
    cH = z(ind.cH(1):ind.cH(2));
        
    % intermediate
    % [cm / bar / s] * ([bar] - [cm^3 bar / micromol / K] * [K] * [micromol / cm^3])
    % = [cm / s]
    Jw = Lp*(delP - (cF - cH).' *(ni.*sigma)*R*T);

    % intermediate
    % [cm / s] * [micromol / cm^2]
    % = [micromol / cm^2 / s]
    Js = B*(cF - cH); 
    
    % consider concentration polarization
    if model_stru.concpolar
        cIn = z(ind.cIn(1):ind.cIn(2));
        
        % intermediate
        % [cm / bar / s] * ([bar] - [cm^3 bar / micromol / K] * [K] * [micromol / cm^3])
        % = [cm / s]
        Jw = Lp*(delP - (cIn - cH).' *(ni.*sigma)*R*T);

        % intermediate
        % [cm / s] * [micromol / cm^2]
        % = [micromol / cm^2 / s]
        Js = B*(cIn - cH);
        
        k = 0.23 * v^0.57 .* D.^0.67 / (nu^0.24 * b^0.43);
        rhs(ind.cIn(1):ind.cIn(2)) = (cF - cH) .* exp(Jw./ k) + cH -cIn;  
    end    

    if data_stru.mode == 'F'
        % dcF_dt
        % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^3] * [cm / s] - [micromol / cm^2 / s])
        % = [micromol / cm^3 / s] 
        rhs(ind.cF(1):ind.cF(2)) = Am * rho / mF * (cF * Jw - Js);
        
    elseif data_stru.mode == 'D'
        % Dialysate concentration
        cD = data_stru.data_config.C_D;
        % dcF_dt
        % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^3] * [cm / s] - [micromol / cm^2 / s])
        % = [micromol / cm^3 / s] 
        rhs(ind.cF(1):ind.cF(2)) = Am * rho / mF * (cD * Jw - Js);

    end  
    
    % dcH_dt, micromol / cm^3 / s  
    % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^2 / s] - [micromol / cm^3] * [cm / s]) 
    % = [micromol / cm^3 / s]
    rhs(ind.cH(1):ind.cH(2)) = Am * rho / mH *(Js - Jw * cH);

    if model_stru.combine_cvmv
        % dmVcV_dt = JwAmrho * ch
        % [micromol / cm^2 / s] * [cm^2] * [g / cm^3]
        % = [micromol * g / cm^3 / s] 
        rhs(ind.cV(1):ind.cV(2)) = Jw * Am * rho * cH;
    else
        % dcV_dt
        % [cm^2] * [g / cm^3] /[g] * ([micromol / cm^2 / s] - [micromol / cm^3] * [cm / s]) 
        % = [micromol / cm^3 / s]
        rhs(ind.cV(1):ind.cV(2)) = - Jw * Am * rho * (- cH + cV) / mV;
    end

end

% dm_dt
% [cm / s] * [cm^2] * [g / cm^3] 
% = [g / s]
dm_dt = - Jw * Am * rho;

% filtration mode
if data_stru.mode == 'F'
    % dmF_dt 
    rhs(ind.mF(1):ind.mF(2)) = dm_dt;
end
% Constant mF for diafiltration mode

% dmV_dt = - dm_dt
% [cm / s] * [cm^2] * [g / cm^3] 
% = [g / s]
rhs(ind.mV(1):ind.mV(2)) = - dm_dt ;

end
