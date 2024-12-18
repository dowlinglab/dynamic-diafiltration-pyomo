function plot_grid_contour(data_stru, model_stru, contour3d_stru)
% Plot grid contours of inividual objectives from 3d search results
% 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - model configuration
%     contour3d_stru - structure save all contours data

x_var = contour3d_stru.var{1};
y_var = contour3d_stru.var{2};
z_var = contour3d_stru.var{3};
xlabelstr = loadlabel(x_var);
ylabelstr = loadlabel(y_var);
zlabelstr = loadlabel(z_var);

grid_density = contour3d_stru.grid_density;

% reshape x,y for contour plot
xx = contour3d_stru.xx;
yy = contour3d_stru.yy;
zz = contour3d_stru.zz;
X = reshape(xx,grid_density,grid_density,6);
Y = reshape(yy,grid_density,grid_density,6);
Z = reshape(zz,grid_density,grid_density,6);

% font size
fs_label = 12;
fs_title = 10;
fs_ax = 10;

% loop over components
for i=1:data_stru.data_config.nc   
    
    
    % reshape objectives for contour plot
    f_m = contour3d_stru.Objectives(i).mass;
    F_m = reshape(f_m,grid_density,grid_density,6);
    
    f_cp = contour3d_stru.Objectives(i).perm_conc;
    F_cp = reshape(f_cp,grid_density,grid_density,6);
    
    if isfield(data_stru.data_raw,'cF_exp')
        f_cr = contour3d_stru.Objectives(i).reten_conc;
        F_cr = reshape(f_cr,grid_density,grid_density,6);
    end
    
    fig = figure;
    for j = 1:size(X,3)  
        % Three plots for individual objectives
        % plot for separated objective - mass residual
        subplot(size(X,3),3,3*(j-1)+1)
        % contour plot
        [c1,h1]=contour(X(:,:,j),Y(:,:,j),F_m(:,:,j),8,'LineWidth',2);
        h1.LevelList = round(h1.LevelList(1:2:end),1);  %rounds levels to 2rd decimal place
        clabel(c1,h1)
        hold     
        set(gca,'FontSize',fs_ax)

        % plot for separated objective - permeate concentration residual
        subplot(size(X,3),3,3*(j-1)+2)
        % contour plot
        [c2,h2]=contour(X(:,:,j),Y(:,:,j),F_cp(:,:,j),'LineWidth',2);
        h2.LevelList = round(h2.LevelList,1);  %rounds levels to 2rd decimal place
        clabel(c2,h2)
        hold
        set(gca,'FontSize',fs_ax)

        % plot for separated objective - retentate concentration residual
        subplot(size(X,3),3,3*(j-1)+3)
        % contour plot
        [c3,h3]=contour(X(:,:,j),Y(:,:,j),F_cr(:,:,j),5,'LineWidth',2);
        h3.LevelList = round(h3.LevelList,1);  %rounds levels to 2rd decimal place
        clabel(c3,h3)        
        hold
        set(gca,'FontSize',fs_ax)
        
        pos = get(subplot(size(X,3),3,3*(j-1)+1),'position');
        annotation('textbox',[0,pos(2),0,0],'string',[zlabelstr,num2str(Z(1,1,j))],'FontSize',fs_label)
    end
    
    % set axes labels and title
    subplot(size(X,3),3,1)
    title({'Log_{10} transform','Mass','Residual Squared[g^2]'},'FontSize',fs_title)
    subplot(size(X,3),3,2)
    title({'Log_{10} transform','Permeate Concentration','Residual Squared[mM^2]'},'FontSize',fs_title)
    subplot(size(X,3),3,3)
    title({'Log_{10} transform','Retentate Concentration','Residual Squared[mM^2]'},'FontSize',fs_title) 
    
    insert_position = regexp(model_stru.titlestr,'filtration');
    modeltitle = [model_stru.titlestr(1:insert_position+10) newline model_stru.titlestr(insert_position+11:end)];
    fi = axes(fig,'visible','off');
    fi.Title.Visible = 'on';
    fi.XLabel.Visible = 'on';
    fi.YLabel.Visible = 'on';
    xlabel(xlabelstr,'FontSize',fs_label)
    ylabel(ylabelstr,'FontSize',fs_label)
    sgtitle({'\fontsize{12}',modeltitle},'interpreter','tex');   
    
    % save figure
    set(gcf,'Units','Inches','Position',[0 0 12 30],'PaperPositionMode','auto')
    print(gcf,['grid_contour-',model_stru.filenamestr,'.pdf'],'-dpdf','-fillpage')

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
            labelstr ='\sigma=';
    end
end
