function plot_sim_only(data_stru, model_stru, sim_stru)
% Plot simulation results only
% 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - partial model configuration
%     sim_stru - simulation results
% 
% Returns:
%     Mass data/prediction comparison plot
%     Concentration data/prediction comparison plot

nc = data_stru.data_config.nc;
figure
% plot mass prediction
subplot(1,nc+1,1)
hold on
% loop over vials
for i = 1:data_stru.data_config.n
    plot(sim_stru(i).time, sim_stru(i).mV,'b-','LineWidth',2);
end
xlabel('Time [s]','FontSize',15)
ylabel('Collected Permeate Mass [g]','FontSize',15);
legend({'Predictions'},'FontSize',12,'Location','Best');
title({'\fontsize{15} Mass Predictions',...
    '\fontsize{11}' model_stru.titlestr},'interpreter','tex');
set(gca,'FontSize',12)
hold off

% plot concentration prediction
% loop over components
for j = 1:nc 
    subplot(1,nc+1,j+1)
    hold on
    % loop over vials
    for i = 1:data_stru.data_config.n
        % retentate concentration predictions
        h1 = plot(sim_stru(i).time, sim_stru(i).cF(:,j),'g-','LineWidth',2);

        % permeate cooncentration predictions
        h2 = plot(sim_stru(i).time(2:end), sim_stru(i).cV(2:end,j),'k-.','LineWidth',2);
        
    end
    legend([h1,h2],{'Retentate (Prediction)','Permeate (Prediction)'},...
        'FontSize',12,'Location','Best');

    xlabel('Time [s]','FontSize',15)
    ylabel(['Concentration of ',data_stru.data_config.namec(j,:),' [mmol/L]'],'FontSize',15);
    title({'\fontsize{15} Concentration Predictions',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');    
    set(gca,'FontSize',12)    
    hold off

end

set(gcf,'Units','Inches','Position',[0.01 0.01 8*(nc+1) 5])
saveas(gcf,['sim_only',data_stru.data_config.namec(j,:),'-',model_stru.filenamestr,'.png'])
end

