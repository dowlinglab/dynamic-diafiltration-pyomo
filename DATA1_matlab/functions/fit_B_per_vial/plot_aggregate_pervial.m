function plot_aggregate_pervial(data_stru, fit_stru)
% Visualize the aggregate B dependence on retentate concentration from multiple experiments

figure
hold on
% plot B per vial vs retentate concentration
for i = 1:length(fit_stru)
    for j = 1:data_stru(i).data_config.n
        cf(j) = fit_stru(i).sim_stru(j).cF(end);
        B(j) = fit_stru(i).B{j};
    end    
    h(i) = plot(cf, B,...
             's','MarkerSize',10);
    set(h(i), 'markerfacecolor', get(h(i), 'color'))
    
    % Delete the paths, leaving only file names. 
    match = "../data_library/";
    exp_name{i} = erase(data_stru(i).filename,match);   
    exp_name{i} = replace(exp_name{i},'_',' ');
end

xlabel('Retentate Concentration [mmol/L]','FontSize',15)
ylabel('B [micrometers/s]','FontSize',15);
title({'\fontsize{15} Salt Permeability vs Retentate Concentration'});

legend(h,exp_name,'Location','eastoutside'); 
set(gcf,'Units','Inches','Position',[0.01 0.01 14 6])
saveas(gcf,['B vs cf_aggregated_',data_stru(1).data_config.namec,'.png'])
end