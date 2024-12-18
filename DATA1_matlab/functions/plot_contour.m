function plot_contour(data_stru, model_stru, contour_stru)
% Plot contours of inividual objectives
% 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - model configuration
%     contour_stru - structure save all contours data

x_var = contour_stru.x_var;
y_var = contour_stru.y_var;
xlabelstr = loadlabel(x_var);
ylabelstr = loadlabel(y_var);

grid_density = contour_stru.grid_density;

% reshape x,y for contour plot
xx = contour_stru.xx;
yy = contour_stru.yy;
X = reshape(xx,grid_density,grid_density);
Y = reshape(yy,grid_density,grid_density);

% font size for labels
fs_label = 15;
fs_title = 15;
fs_ax = 12;

% loop over components
for i=1:data_stru.data_config.nc   
    
    % reshape objectives for contour plot
    f_m = contour_stru.Objectives(i).mass;
    ind_f_m = contour_stru.Objectives(i).m_min_ind;
    F_m = reshape(f_m,grid_density,grid_density);
    
    f_cp = contour_stru.Objectives(i).perm_conc;
    ind_f_cp = contour_stru.Objectives(i).cp_min_ind;
    F_cp = reshape(f_cp,grid_density,grid_density);
    
    if isfield(data_stru.data_raw,'cF_exp')
        f_cr = contour_stru.Objectives(i).reten_conc;
        ind_f_cr = contour_stru.Objectives(i).cr_min_ind;
        F_cr = reshape(f_cr,grid_density,grid_density);
    end
    
    fig = figure;
    
    % Three plots for individual objectives
    % plot for separated objective - mass residual
    subplot(1,3,1)
    % contour plot
    contour(X,Y,F_m,'ShowText','on','LineWidth',2);
    hold
    % add minimum point to plot
    mini_pt_m = plot(xx(ind_f_m),yy(ind_f_m),'^','markersize', 8,...
        'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
    legend(mini_pt_m,['Minimum(',num2str(xx(ind_f_m)), ',',num2str(yy(ind_f_m)), ',',num2str(min(f_m)), ')'],...
        'Location','southoutside'); 
    ylabel(ylabelstr,'FontSize',fs_label)
    set(gca,'FontSize',fs_ax)
    title({'Log_{10} transform','Mass','Residual Squared[g^2]'},'FontSize',fs_title)

    % plot for separated objective - permeate concentration residual
    subplot(1,3,2)
    % contour plot
    contour(X,Y,F_cp,'ShowText','on','LineWidth',2);
    hold
    % add minimum point to plot
    mini_pt_cp = plot(xx(ind_f_cp),yy(ind_f_cp),'^','markersize', 8,...
        'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
    legend(mini_pt_cp,['Minimum(',num2str(xx(ind_f_cp)), ',',num2str(yy(ind_f_cp)), ',',num2str(min(f_cp)), ')'],...
        'Location','southoutside');
    set(gca,'FontSize',fs_ax)
    title({'Log_{10} transform','Permeate Concentration','Residual Squared[mM^2]'},'FontSize',fs_title)

    % plot for separated objective - retentate concentration residual
    subplot(1,3,3)
    % contour plot
    contour(X,Y,F_cr,'ShowText','on','LineWidth',2);
    hold
    % add minimum point to plot
    mini_pt_cr = plot(xx(ind_f_cr),yy(ind_f_cr),'^','markersize', 8,...
        'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
    legend(mini_pt_cr,['Minimum(',num2str(xx(ind_f_cr)), ',',num2str(yy(ind_f_cr)), ',',num2str(min(f_cr)), ')'],...
        'Location','southoutside');
    set(gca,'FontSize',fs_ax)
    title({'Log_{10} transform','Retentate Concentration','Residual Squared[mM^2]'},'FontSize',fs_title)

    % set axes labels and title
    xlabel(axes(fig,'visible','off'),xlabelstr,'FontSize',fs_label,'visible','on')
    sgtitle({'\fontsize{18} Log Transform Squared Norm of Residuals',...
        '\fontsize{12}',model_stru.titlestr},'interpreter','tex');   

    % save figure
    set(gcf,'Units','Inches','Position',[2 2 14 7],'Papersize',[60 20],'PaperPositionMode','auto')
    saveas(gcf,[model_stru.filenamestr,'\objcontour-','x_',x_var,'-y_',y_var,'.png'])

end
end

function labelstr = loadlabel(var_name)
%load axis label name
    switch var_name
        case "Lp"
            labelstr = 'L_p   [L / m^{2} / h / bar]';
        case "B"
            labelstr = 'B   [\mum / s]';
        case "sigma"
            labelstr ='\sigma   [dimensionless]';
    end
end
