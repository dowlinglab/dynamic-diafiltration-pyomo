function plot_pervial_vsclassical(data_stru, model_stru, fit_stru)
% Plot B dependence on retentate concentration comparing with classical
% analysis

data_pre = classical_analysis(data_stru.dataset,false,true,true);

% plot B per vial vs retentate concentration
% loop over components
for j = 1:data_stru.data_config.nc 
    figure
    hold on
    for i = 1:data_stru.data_config.n
        cf(i) = fit_stru.sim_stru(i).cF(end,j);
        B(i) = fit_stru.B{i}(j,j);
    end    
    h1 = plot(cf, B,...
             'cs','MarkerSize',10,'DisplayName','Regression');
    
    h2 = plot(data_pre.cf(:,j), data_pre.B(:,j),...
         's','MarkerSize',10,'DisplayName','Classical analysis');
    set(h2, 'markerfacecolor', get(h2, 'color'))    
         
    xlabel(['Retentate concentration of ',data_stru.data_config.namec(j,:),' [mmol/L]'],'FontSize',15);
    ylabel('B [micrometers/s]','FontSize',15)
    title({'\fontsize{15} Salt Permeability vs Retentate Concentration',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');
    legend
    set(gca,'FontSize',12)

    hold off
    saveas(gcf,[model_stru.filenamestr,'\B vs cf comparison ',data_stru.data_config.namec(j,:),'-',model_stru.filenamestr,'.png'])
end
end