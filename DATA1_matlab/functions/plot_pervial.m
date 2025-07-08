function plot_pervial(data_stru, model_stru, fit_stru)
% Plot B dependence on retentate concentration

% plot B per vial vs retentate concentration
% loop over components
for j = 1:data_stru.data_config.nc 
    figure
    hold on
    for i = 1:data_stru.data_config.n
        h1 = plot(fit_stru.sim_stru(i).cF(end,j),fit_stru.B{i}(j,j),'cs','MarkerFaceColor','c','MarkerSize',10);
    end
    xlabel(['Retentate concentration of ',data_stru.data_config.namec(j,:),' [mmol/L]'],'FontSize',15);
    ylabel('B [micrometers/s]','FontSize',15)
    title({'\fontsize{15} Salt Permeability vs Retentate Concentration',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');
    set(gca,'FontSize',12)

    hold off
    saveas(gcf,[model_stru.filenamestr,'\B vs cf ',data_stru.data_config.namec(j,:),'-',model_stru.filenamestr,'.png'])
end
end