function sim_stru = sim_model(theta, data_stru, model_stru)
% Save simulation results for each vial into a structure
% 
% Arguments:   
%     theta - fitting parameters
%     data_stru - all experimental data
%     model_stru - partial model configuration
% 
% Returns:
%     sim_stru - simulation results which includes
%         time - simulation time span
%         cF - retentate concentration
%         mF - retentate mass
%         cV - permeate concentration
%         mV - permeate mass
%         cH - holdup concentration

% variable indices
ind = model_stru.var_ind;

% Integrate variable initial values
M_V0 = model_stru.initialization.M_V0;

% Toggle for filtration/diafiltration model 
if data_stru.mode == 'F'
    z0 = [model_stru.initialization.M_F0;model_stru.initialization.C_F0;M_V0];
elseif data_stru.mode == 'D'
    % constant mF therefore not an integrated variable in diafiltration
    z0 = [model_stru.initialization.C_F0;M_V0];
end    

if model_stru.combine_cvmv
    C_Vint0 = model_stru.initialization.C_V0 * M_V0;
else
    C_Vint0 = model_stru.initialization.C_V0;
end
z0 = [z0; C_Vint0];

% Use theta to update model structure
model_stru = update_model_stru(theta, model_stru);

if model_stru.holdup
    % if initial hold up concentration chi is a fitting parameter(check for all components)
    if ~all(model_stru.ch0_fixed)
        z0 = [z0; model_stru.model_para.chi];  
    % else all chi values are fixed
    else
        z0 = [z0; model_stru.initialization.C_H0];
    end
end

my_ode = @(t,z) diafiltration_ode(t, z, data_stru, model_stru);

if model_stru.concpolar   
    z0 = [z0; model_stru.initialization.C_F0];
    M = eye(length(z0));
    M(ind.cIn(1):ind.cIn(2),ind.cIn(1):ind.cIn(2)) = 0;
else
    M = eye(length(z0));
end
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-10);

% Loop over vials
for i = 1:data_stru.data_config.n
    
    tspan = data_stru.data_raw(i).time;
    % Simulate
    [T, Z1] = ode15s(my_ode,tspan,z0,opts);
    
    % store results
    sim_stru(i).time = T;
    
    if data_stru.mode == 'F'
        sim_stru(i).mF = Z1(:,ind.mF(1):ind.mF(2)); % g
    end    
    
    sim_stru(i).cF = Z1(:,ind.cF(1):ind.cF(2)); % mmol/L = micromol / cm^3
    
    % Reset all vial mass with zero initial value
    sim_stru(i).mV = Z1(:,ind.mV(1):ind.mV(2)) - M_V0; % g
    
    if model_stru.combine_cvmv
        % the product cv*mv is integrated as a single state and returned as
        % a product by the ode function. cv(t) = cv*mv(t)/mv(t), t>0
        sim_stru(i).cV = Z1(:,ind.cV(1):ind.cV(2))./Z1(:,ind.mV(1):ind.mV(2)); % mmol/L
    else
        % cv is integrated as an independent state and returned directly
        % from the ode function
        sim_stru(i).cV = Z1(:,ind.cV(1):ind.cV(2)); % mmol/L
    end
    
    if model_stru.holdup
        sim_stru(i).cH = Z1(:,ind.cH(1):ind.cH(2)); % mmol/L
    end
    
    % advance initial condition
    C_F0 = Z1(end,ind.cF(1):ind.cF(2)).';
    C_V0 = sim_stru(i).cV(end,:).';
    
    if model_stru.holdup
        C_H0 = Z1(end,ind.cH(1):ind.cH(2)).';
        C_V0 = C_H0; % initialize vial conc. by holdup conc.
    end
    
    if model_stru.combine_cvmv
        C_Vint0 = C_V0 * M_V0;
    else
        C_Vint0 = C_V0;
    end 
    
    if model_stru.concpolar    
        C_In0 = Z1(end,ind.cIn(1):ind.cIn(2)).';
    end
    
    % assemble all initial conditions
    if data_stru.mode == 'F'
        M_F0 = Z1(end,ind.mF(1):ind.mF(2)).';
        z0 = [M_F0; C_F0; M_V0; C_Vint0];
    elseif data_stru.mode == 'D'
        % constant mF in diafiltration
        z0 = [C_F0; M_V0; C_Vint0];
    end  
   
    if model_stru.holdup
        z0 = [z0; C_H0];
    end
    
    if model_stru.concpolar    
        z0 = [z0; C_In0];
    end
end

% save(['sim_stru-Dataset', num2str(data_stru.dataset), 'model', num2str(model_stru.filenamestr), '.mat'], 'sim_stru')
end

