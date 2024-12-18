function plot_sim(data_stru, model_stru, sim_stru)
% Plot simulation results comparing with data
% 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - partial model configuration
%     sim_stru - simulation results
% 
% Returns:
%     Mass data/prediction comparison plot
%     Concentration data/prediction comparison plot

% plot mass data/prediction comparison 
figure
hold on
t_delay = data_stru.data_raw(1).time(1);
% loop over vials
for i = 1:data_stru.data_config.n
    plot(data_stru.data_raw(i).time-t_delay, data_stru.data_raw(i).mass,'r.');
    plot(sim_stru(i).time-t_delay, sim_stru(i).mV,'b-','LineWidth',2);
end
xlabel('Time [s]','FontSize',15)
ylabel('Collected Permeate Mass [g]','FontSize',15);
legend({'Measurements','Predictions'},'FontSize',12,'Location','Best');
title({'\fontsize{15} Mass Predictions',...
    '\fontsize{11}' model_stru.titlestr},'interpreter','tex');
set(gca,'FontSize',12)
hold off

saveas(gcf,[model_stru.filenamestr,'\mass-',model_stru.filenamestr,'.png'])

% plot concentration data/prediction comparison
% loop over components
for j = 1:data_stru.data_config.nc 
    figure
    hold on
    
    if isfield(data_stru.data_raw(end),'cF_exp')
        % retentate concentration measurements
        for i = 1:data_stru.data_config.n
            % if retentate concentration is measured for each vial
            if ~isempty(data_stru.data_raw(i).cF_exp)
                % one measurement for each vial
                if size(data_stru.data_raw(i).cF_exp,1) == 1
                    plot(data_stru.data_raw(i).time(end)-t_delay,data_stru.data_raw(i).cF_exp(j),...
                        'ms','MarkerFaceColor','m','MarkerSize',10);
                % continous measurements
                else
                    plot(data_stru.data_raw(i).time-t_delay, data_stru.data_raw(i).cF_exp(:,j),...
                        'mo','MarkerFaceColor','w','MarkerSize',5);
                end
            end
        end     
    end
    
    % loop over vials
    for i = 1:data_stru.data_config.n
        % retentate concentration predictions
        h1 = plot(sim_stru(i).time-t_delay, sim_stru(i).cF(:,j),'g-','LineWidth',2);

        % vial concentration predictions
        h2 = scatter(sim_stru(i).time(end)-t_delay, sim_stru(i).cV(end,j),100,...
            'ks','MarkerFaceColor','k','MarkerFaceAlpha',0.2);
        
        if isfield(sim_stru,'cH')
            h2h = plot(sim_stru(i).time(:)-t_delay, sim_stru(i).cH(:,j),'k-','LineWidth',2);
        end
        
        % vial concentration measurements
        h3 = scatter(data_stru.data_raw(i).time(end)-t_delay, data_stru.data_raw(i).cV_avg(j),100,...
            'cs','MarkerFaceColor','c','MarkerFaceAlpha',0.2);
    end
    
    % initial feed concentration before overflow
    h4 = plot(0 - t_delay,data_stru.data_config.C_F0(j),...
        'ms','MarkerFaceColor','m','MarkerSize',10);
    % initial feed concentration after overflow
    plot(0,model_stru.initialization.C_F0(j),...
        'ms','MarkerFaceColor','m','MarkerSize',10);
    
%     if isfield(sim_stru,'cH')
%         legend([h1,h4,h2h,h2,h3],{'Retentate (Prediction)','Retentate (Measurement)','Permeate (Prediction)','Vial (Prediction)','Vial (Measurement)'},...
%             'FontSize',12,'Location','Best');
%     else
%         legend([h1,h4,h2,h3],{'Retentate (Prediction)','Retentate (Measurement)','Vial (Prediction)','Vial (Measurement)'},...
%             'FontSize',12,'Location','Best');
%     end

    xlabel('Time [s]','FontSize',15)
    ylabel(['Concentration of ',data_stru.data_config.namec(j,:),' [mM]'],'FontSize',15);
    title({'\fontsize{15} Concentration Predictions',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');    
    set(gca,'FontSize',12)
    
    hold off
    saveas(gcf,[model_stru.filenamestr,'\concentration ',data_stru.data_config.namec(j,:),'-',model_stru.filenamestr,'.png'])
end

end

