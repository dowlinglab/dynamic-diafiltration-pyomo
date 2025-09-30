function plot_sim_v2(data_stru, model_stru, sim_stru, opts)

% Function to plot simulation results against experimental data

% Inputs:
%   - data_stru - structured experimental data
%   - model_stru - structure with partial model configuration (e.g., title, file names)
%   - sim_stru - structure with simulation results (time, concentration, mass)
%   - opts - options structure controlling plot behavior

% Returns:
%     Mass data/prediction comparison plot
%     Concentration data/prediction comparison plot

% Plot mass data/prediction comparison

% If opts is not provided, initialize with empty struct; incase opts is not
% user-defined
if nargin < 4
    opts = struct();
end

% Set default behavior for saving, showing, and layout
if ~isfield(opts, 'do_save'), opts.do_save = true; end              % Save plots to file?
if ~isfield(opts, 'do_show'), opts.do_show = true; end              % Show figures?
if ~isfield(opts, 'subplot_mode'), opts.subplot_mode = false; end   % Use subplot layout?
if ~isfield(opts, 'fig_id_mass'), opts.fig_id_mass = figure; end    % Mass plot figure handle

% Compute delay from experiment start
t_delay = data_stru.data_raw(1).time(1);
n_vials = data_stru.data_config.n;     % Number of vials
n_comp = data_stru.data_config.nc;     % Number of chemical components

% === MASS PLOT ===
if opts.do_show
    figure(opts.fig_id_mass); % Activate the specified figure for mass plot
end
hold on

% Plot experimental and simulated mass for each vial
for i = 1:n_vials
    plot(data_stru.data_raw(i).time - t_delay, data_stru.data_raw(i).mass, 'r.');     % Experimental data
    plot(sim_stru(i).time - t_delay, sim_stru(i).mV, 'b-', 'LineWidth', 2);           % Simulated data
end

% Label and format mass plot
xlabel('Time [s]', 'FontSize', 15)
ylabel('Collected Permeate Mass [g]', 'FontSize', 15)
legend({'Measurements', 'Predictions'}, 'FontSize', 12, 'Location', 'Best')
title({'\fontsize{15} Mass Predictions', ['\fontsize{11}' model_stru.titlestr]}, 'interpreter', 'tex')
set(gca, 'FontSize', 12)
hold off

% Save figure if enabled
if opts.do_save
    saveas(gcf, fullfile(model_stru.filenamestr, ['mass-', model_stru.filenamestr, '.png']));
end

% === CONCENTRATION PLOTS ===
% plot concentration data/prediction comparison
% loop over components

for j = 1:n_comp
    % Determine layout for concentration plot
    if opts.subplot_mode
        if j == 1
            fig_conc = figure; % Create new figure for subplot
        end
        subplot(n_comp, 1, j); % Create subplot for each component
    elseif isfield(opts, 'fig_id_conc')
        figure(opts.fig_id_conc(j)); % Use external figure handle
    else
        figure; % Create new figure per component
    end

    hold on

    % Plot retentate concentration measurements (if available)
    if isfield(data_stru.data_raw(end), 'cF_exp')
        for i = 1:n_vials % retentate concentration measurements
            if ~isempty(data_stru.data_raw(i).cF_exp) % if retentate concentration is measured for each vial
                if size(data_stru.data_raw(i).cF_exp, 1) == 1 % one measurement for each vial
                    % One data point
                    plot(data_stru.data_raw(i).time(end) - t_delay, data_stru.data_raw(i).cF_exp(j), ...
                        'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 10); % continous measurements
                else
                    % Time-series data
                    plot(data_stru.data_raw(i).time - t_delay, data_stru.data_raw(i).cF_exp(:, j), ...
                        'mo', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
                end
            end
        end
    end

    % Plot predicted concentrations and vial data
    for i = 1:n_vials
        plot(sim_stru(i).time - t_delay, sim_stru(i).cF(:, j), 'g-', 'LineWidth', 2); % Retentate conc. prediction
        scatter(sim_stru(i).time(end) - t_delay, sim_stru(i).cV(end, j), 100, ...      % Vial concentration prediction
            'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2);

        if isfield(sim_stru, 'cH')
            plot(sim_stru(i).time - t_delay, sim_stru(i).cH(:, j), 'k-', 'LineWidth', 2); % Optional: permeate concentration
        end

        scatter(data_stru.data_raw(i).time(end) - t_delay, data_stru.data_raw(i).cV_avg(j), 100, ... % Measured vial concentration
            'cs', 'MarkerFaceColor', 'c', 'MarkerFaceAlpha', 0.2);
    end

    % Initial concentrations before and after overflow
    plot(0 - t_delay, data_stru.data_config.C_F0(j), 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
    plot(0, model_stru.initialization.C_F0(j), 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 10);

    % Axis labels and title
    xlabel('Time [s]', 'FontSize', 15)
    ylabel(['Concentration of ', strtrim(data_stru.data_config.namec(j,:)), ' [mM]'], 'FontSize', 15)
    title({'\fontsize{15} Concentration Predictions', ['\fontsize{11}' model_stru.titlestr]}, 'interpreter', 'tex')
    set(gca, 'FontSize', 12)
    hold off

    % Save concentration plot if enabled
    if opts.do_save && ~opts.subplot_mode
        saveas(gcf, fullfile(model_stru.filenamestr, ...
            ['concentration_', strtrim(data_stru.data_config.namec(j,:)), '-', model_stru.filenamestr, '.png']));
    end
end
end

