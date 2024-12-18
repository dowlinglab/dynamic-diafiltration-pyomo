function data_stru = load_data(dataset_number)
% Load all the data into a structure data_stru
% 
% Arguments:
%     dataset_number - specify data set e.g., 1.1, 1.2, etc.
% 
% Returns:
%     data_stru - a structrue which includes
%         dataset - dataset selection [float]
%         filename - filename prefix corresponds to dataset in data_library
%         mode - experiement mode
%         continuous_cF - continuous retentate concentration measurement toggle
%         data_config - variables & parameters initialization      
%             M_F0 - initial feed mass [g]
%             C_F0 - initial feed concentration [mMol/L]
%             M_O - overflow mass [g]
%             C_D - dialysate concentration [mMol/L]
%             nc - number of components [-]
%             namec - names of components [-]
%             ni - number of dissolved species [-]
%             delP - applied pressure [bar]
%             Temp - temperature [K]
%             Am - effective membrane area [cm^2]
%             rho - density [g/cm^3]
%             Lp0 - initial guess for Lp [L/m^2/bar/h]
%             B0 - initial guess for B [um/s]
%             sigma0 - initial guess for sigma [-]
%             theta0 - assembles all initial guess
%             n - total number of vials [-] 
%             nr - total number of residuals [-]
%         data_raw - Experimental raw data, contains array structure of each vial
%             number - sorting number of vials [-]
%             time - time data [s]
%             mass - mass measurements [g]
%             cV_avg - permeate concentration measurements [mMol/L]
%             cF_exp - final retentate concentration measurement [mMol/L]
%             nr - number of residuals [-]

data_stru.dataset = dataset_number;

% Read in data dictionary to map shorthand code(e.g., dataset) to
% filename and experiment mode
load_data_dict(dataset_number)
% Read in experiment configuration
load_configuration_file(data_stru.filename)
% Read in experiment measurements
load_raw_data(data_stru.filename)

% number of objectives
% mass and permeate concentration objectives
data_stru.n_obj = size(data_stru.data_config.M_F0,2) ...
                + size(data_stru.data_config.C_F0,2);
% retentate concentration objectives       
if isfield(data_stru.data_raw,'cF_exp')
    data_stru.n_obj = data_stru.n_obj + size(data_stru.data_config.C_F0,2); 
end

% These functions can directly modify dat_struct
    function load_data_dict(dataset_number)
        % Read in data dictionary to map shorthand code(e.g., dataset) to
        % filename and experiment mode
        T = readtable('../data_library/data_dictionary.csv');
        filename = char(T.Filename(T.Index == dataset_number));
        data_stru.filename = ['../data_library/',filename];
        data_stru.mode = char(T.Mode(T.Index == dataset_number));
    end

    function load_configuration_file(filename)
        % Experiment configuration
        T = readtable([filename,'_config.csv']);
        % Specify initial conditions 
        % Initial mass of feed (initial filtation volume is 10ml)        
        data_stru.data_config.M_F0 = T.initialFeedMass(~isnan(T.initialFeedMass)); %[g] 
        % Initial concentration of feed
        data_stru.data_config.C_F0 = T.initialFeedConcentration; % [mmol/l]
        
        if data_stru.mode == 'D'
            % Overflow mass 
            data_stru.data_config.M_O = T.overflowMass; %[g] 
            % Dialysate concentration
            data_stru.data_config.C_D = T.dialysateConcentration; % [mmol/l]
        end
        
        % Number of components
        nc = length(data_stru.data_config.C_F0);
        data_stru.data_config.nc = nc;
        data_stru.data_config.namec = char(T.component);
        % Number of dissolved species
        data_stru.data_config.ni = T.numberOfDissolvedSpecies;
        % Applied pressure
        data_stru.data_config.delP = T.appliedPressure(~isnan(T.appliedPressure));
        % Temperature
        data_stru.data_config.Temp = T.temperature(~isnan(T.temperature));
        % Effective membrane area
        data_stru.data_config.Am = T.effectiveMembraneArea(~isnan(T.effectiveMembraneArea));
        % Density 
        data_stru.data_config.rho = T.density(~isnan(T.density));     
        
        % Initial guess of parameters
        data_stru.data_config.Lp0 = T.initialLp(~isnan(T.initialLp)); 
        Bidx = find(string(T.Properties.VariableNames) == "initialB");
        data_stru.data_config.B0 = table2array(T(:,Bidx:(Bidx + nc -1)));
        data_stru.data_config.sigma0 = T.initialSigma; 
        
        % parameter order: Lp, B(in column), sigma
        data_stru.data_config.theta0 = [data_stru.data_config.Lp0;... 
                                        reshape(data_stru.data_config.B0,[],1);...
                                        data_stru.data_config.sigma0];
    end

    function load_raw_data(filename)
        % Read in experiment measurements
        exp_index = readmatrix([filename,'_ind.csv']);
        time_table = readtable([filename,'_m.csv']);
        time = time_table{:,1};
        mass = time_table{:,2};
        cV_exp = readmatrix([filename,'_cV.csv']);
        cF_table = readtable([filename,'_cF.csv']);
        cF = cF_table{:,:};
        % continuous retentate concentration measurement toggle
        data_stru.conductivity_cF = true;
   
        % number of experiments
<<<<<<< HEAD
        data_stru.data_config.n = size(exp_index,1);
        % number of residuals
        data_stru.data_config.nr = 0;
 
=======
        n = size(exp_index,1);
        n_extra = size(exp_index,1)-exp_index(end,1);
        data_stru.data_config.n = n;
        data_stru.data_config.n_extra = n_extra;
        data_stru.data_config.n_v0 = 1;
        if n_extra ~= 0
            update_exp_index = exp_index(1:n_extra,:);
            [sd,r] = sort(exp_index(1:n_extra,3),'ascend');
            update_exp_index(2:3,2) = sd(1:2)';
            update_exp_index(1:3,3) = sd(1:3)';
            exp_index(1:3,:) = update_exp_index;
            if exp_index(3,3) ~= exp_index(4,2)
                exp_index(4,2) = sd(3);
                data_stru.data_config.n_v0 = find(r==3)+1;
            else 
                data_stru.data_config.n_v0 = n_extra+1;
            end
            data_stru.data_config.n_h = find(r==1);
            data_stru.data_config.n_A = find(r==2);
        end

        % number of residuals
        data_stru.data_config.nr = 0;
        
        % vial 0 time extrapolation
        if data_stru.data_config.n_extra ~= 0
            A=[];
            t=[];
            for i = data_stru.data_config.n_extra + 1 : n;
                jstart = exp_index(i,2);
                jend = exp_index(i,3);
                ind = jstart:jend;
                A = [A; mass(ind)];
                t = [t; time(ind)];
            end
            A = fillmissing(A,'linear');
            p = polyfit(A,t,6);
            time_vial0 = polyval(p,-0.68)
%             t_ = time(exp_index(1,2)):5:time(exp_index(data_stru.data_config.n,3));
            a = mass(exp_index(data_stru.data_config.n_extra + 1,2))-0.8:0.01:mass(exp_index(data_stru.data_config.n_extra + 1,2));
            A_ = [a';A];
            t_ = polyval(p,A_);
            figure
            plot(t,A,'r')
            hold on
            plot(t_,A_,'b')
%             time_vial0 = interp1(A,time(ind),mass(jstart)-0.8, 'nearest', 'extrap')
%             mass_vial0 = interp1(t,A,200,'spline','extrap')
        end
>>>>>>> lag&overflow
        % loop over vials, store data into a structured array
        for i = 1 : n
            % store vial number
            data_stru.data_raw(i).number = i - data_stru.data_config.n_extra;
 
            % grab time index
            jstart = exp_index(i,2);
            jend = exp_index(i,3);
            ind = jstart:jend;

            % store time - lag time removed
            data_stru.data_raw(i).time = time(ind);

            % store mass of vial
            data_stru.data_raw(i).mass = mass(ind) - mass(jstart);

            % eliminate any gross outliers
            data_stru.data_raw(i).mass(data_stru.data_raw(i).mass < 0) = NaN;
            % data_stru.data_raw(i).mass(data_stru.data_raw(i).mass > 10) = NaN;
            data_stru.data_raw(i).mass(data_stru.data_raw(i).mass > 1.2) = NaN;
            
            nr = sum(~isnan(data_stru.data_raw(i).mass));
            
            if cF_table.Properties.VariableNames{1} == "time"
                data_stru.conductivity_cF = true;
                % store continous retentate conc. measurements
                data_stru.data_raw(i).cF_exp = cF(ind,2:end);
                
                % eliminate any gross outliers
                data_stru.data_raw(i).cF_exp(data_stru.data_raw(i).mass < 0) = NaN;
                data_stru.data_raw(i).cF_exp(data_stru.data_raw(i).mass > 1.2) = NaN;
                
                nr = nr + sum(~isnan(data_stru.data_raw(i).cF_exp(:,1)));
            end

<<<<<<< HEAD
            % store vial concentration
            data_stru.data_raw(i).cV_avg = cV_exp(i,2:end);
            % determine number of residuals
            data_stru.data_raw(i).nr = nr + size(cV_exp,2) - 1;           
            data_stru.data_config.nr = data_stru.data_config.nr + data_stru.data_raw(i).nr;

=======
            if cV_table.Properties.VariableNames{1} == "time"
                % corrected time for vial concentration cV
                data_stru.data_raw(i).cV_avg = NaN(length(data_stru.data_raw(i).time),width(cV_table)-1);
                data_stru.data_raw(i).cV_avg(ismember(data_stru.data_raw(i).time, cV(:,1))) = cV(ismember(cV(:,1),data_stru.data_raw(i).time),2:end);
            else
                % store vial concentration
                data_stru.data_raw(i).cV_avg = cV_exp(i,2:end);
                % determine number of residuals
                data_stru.data_raw(i).nr = nr + size(cV_exp,2) - 1;           
                data_stru.data_config.nr = data_stru.data_config.nr + data_stru.data_raw(i).nr;
            end
>>>>>>> lag&overflow
        end

        % retentate conc. measurements
        % if retentate conc. is measured
        if size(cF,1) > 1 && cF_table.Properties.VariableNames{1} == "numberOfVial"
            % only final retentate conc. is measured
            if cF(end,1) == -1
                data_stru.data_raw(end).cF_exp = cF(end,2:end);
                data_stru.data_config.nr = data_stru.data_config.nr + size(cF,2) - 1;
            % retentate conc. measured for each vial
            else
                data_stru.conductivity_cF = true;
                for i = 1:n
                    data_stru.data_raw(i).cF_exp = cF(i+1,2:end);% from 0(initial)
                end
            end 
        end       
    end

% save(['data_stru-dataset', num2str(dataset_number), '.mat'], 'data_stru')
end

