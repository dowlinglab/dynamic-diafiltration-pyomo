% Generate heatmap for sigma sensitivity on different experiment conditions
clear all
Dat = 511.12;
cf0 = 5;
cd = 1:1:80;
End_vial = [1;5;10];
delp = 10:1:120;
% Sigma values evaluated at
sigma = [0.9,1];
% Parameters evaluated at
theta = [3.8994;0.29489;1];% NF90.5 diafiltration
% theta = [];

log_switch = false;
view_heatmap = false;
osmotic_switch = true;
Jw0_filter = false;
Jw_filter = true;
scale_switch = true;

for e = 1:length(End_vial)
end_vial = End_vial(e);

% Run simulation
data_stru = load(['data/data_stru-dataset', num2str(Dat), '.mat']).data_stru; 
legacy_model_name = '201cvmv';

% Initialization
range_mV = zeros(length(delp),length(cd));
range_cF = zeros(length(delp),length(cd));
range_cH = zeros(length(delp),length(cd));
range_pO = zeros(length(delp),length(cd));
pO_values = zeros(length(delp),length(cd),length(sigma));


for c = 1:length(cd)
    C_D = cd(c);        
    data_stru.data_config.C_D = C_D;
    
    for p = 1:length(delp)
        delP = delp(p);
        data_stru.data_config.delP = delP/14.504;

        model_stru = model_config(legacy_model_name,data_stru);
        model_stru.concpolar = true;
        model_stru = update_model_config(data_stru,model_stru);
        model_stru.initialization.C_F0 = cf0;
        
        % Simulate
        sigsen_stru = sigma_sensitivity(data_stru, model_stru,sigma,theta,false);
        mV = zeros(length(sigsen_stru),1);
        cF = zeros(length(sigsen_stru),1);
        cH = zeros(length(sigsen_stru),1);

        nRT = data_stru.data_config.ni * 8.314e-5 * data_stru.data_config.Temp;
        % save ending values for variables
        for s = 1:length(sigsen_stru)
            mV(s) = sigsen_stru(s).sim_stru(end_vial).mV(end);
            cF(s) = sigsen_stru(s).sim_stru(end_vial).cF(end);
            cH(s) = sigsen_stru(s).sim_stru(end_vial).cH(end);
            CF_all=[];
            CH_all=[];
            for v = 1:data_stru.data_config.n
                CF_all = [CF_all;sigsen_stru(s).sim_stru(v).cF];
                CH_all = [CH_all;sigsen_stru(s).sim_stru(v).cH];                
            end
            Osmotic(s).cF = CF_all;
            Osmotic(s).cH = CH_all;
            Osmotic(s).pO = nRT * (CF_all - CH_all);              
        end
        
        % Jw(0)
        if Jw0_filter
            Jw0_Lp = data_stru.data_config.delP - nRT * sigma * C_F0;
            if any(Jw0_Lp<0)
                range_mV(p,c) = nan;
                range_cF(p,c) = nan;
                range_cH(p,c) = nan;
                range_pO(p,c) = nan;
                continue
            end
        end

        % Jw
        if Jw_filter
            for v = 1:data_stru.data_config.n
                for s = 1:length(sigsen_stru)
                    cF_max(s) = max(sigsen_stru(s).sim_stru(v).cF);
                end
                Jw_Lp(v) = data_stru.data_config.delP - nRT * max(sigma) * max(cF_max);
            end
            
            if any(Jw_Lp<1e-8)
                range_mV(p,c) = nan;
                range_cF(p,c) = nan;
                range_cH(p,c) = nan;
                range_pO(p,c) = nan;
                continue
            end
        end
        
        % Osmotic pressure
        pO = nRT * (cF - cH);
        
        for s=1:length(sigma)
            pO_values(p,c,s) = pO(s);
        end
        
        % log transformed
        if log_switch
            range_mV(p,c) = log10(range(mV));
            range_cF(p,c) = log10(range(cF));
            range_cH(p,c) = log10(range(cH));
            range_pO(p,c) = log10(range(pO));
        else
            range_mV(p,c) = range(mV);
            range_cF(p,c) = range(cF);
            range_cH(p,c) = range(cH);
            range_pO(p,c) = range(pO);
        end
        
        % scaled by uncertainty
        if scale_switch
            range_mV(p,c) = range_mV(p,c)/0.01;
            range_cF(p,c) = range_cF(p,c)/(mean(cF)*0.003);
            range_cH(p,c) = range_cH(p,c)/(mean(cH)*0.03);
        end         
    end
end

% save results
contour_sig_stru.cd = cd;
contour_sig_stru.delp = delp;
contour_sig_stru.range_mV = range_mV;
contour_sig_stru.range_cF = range_cF;
contour_sig_stru.range_cH = range_cH;
stru_name = ['contour_sig_stru-dat', num2str(Dat), '-endvial', num2str(end_vial),'.mat'];
save(stru_name, 'contour_sig_stru')
movefile(stru_name,'sigma_sensitivity')

% Visualization
if view_heatmap
    figure(1)
    subplot(1,3,1)
    h1=heatmap(cd,fliplr(delp),flipud(range_mV));
    h1.Title = '\fontsize{15} Difference in Mass Predictions';
    h1.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h1.YLabel = '\fontsize{15} Applied pressure [psi]';
    h1.Colormap = parula;
    h1.ColorbarVisible = 'off';
    h1.MissingDataColor = [0.5, 0.5, 0.5];
    h1.GridVisible = 'off';


    subplot(1,3,2)
    h2=heatmap(cd,fliplr(delp),flipud(range_cF));
    h2.Title = '\fontsize{15} Difference in Retentate Predictions';
    h2.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h2.YLabel = '\fontsize{15} Applied pressure [psi]';
    h2.Colormap = parula;
    h2.ColorbarVisible = 'off';
    h2.MissingDataColor = [0.5, 0.5, 0.5];
    h2.GridVisible = 'off';

    subplot(1,3,3)
    h3=heatmap(cd,fliplr(delp),flipud(range_cH));
    h3.Title = '\fontsize{15} Difference in Permeate Predictions';
    h3.XLabel = '\fontsize{15} Dialysate concentration [mM]';
    h3.YLabel = '\fontsize{15} Applied pressure [psi]';
    h3.Colormap = parula;
    h3.ColorbarVisible = 'off';
    h3.MissingDataColor = [0.5, 0.5, 0.5];
    h3.GridVisible = 'off';

    if length(cd) <  11
        set(gcf,'Units','Inches','Position',[0.01 0.01 24 5])
    else
        set(gcf,'Units','Inches','Position',[0.01 0.01 48 10])
    end

    plotname = ['sigma_sensitivity_heatmap-',model_stru.filenamestr,'.png'];
    saveas(gcf,plotname)
    movefile(plotname,'sigma_sensitivity')

    if osmotic_switch
        figure(2)
        subplot(1,1+length(sigma),1)
        h1=heatmap(cd,fliplr(delp),flipud(range_pO));
        h1.Title = '\fontsize{15} Difference in Osmotic Pressure';
        h1.XLabel = '\fontsize{15} Dialysate concentration [mM]';
        h1.YLabel = '\fontsize{15} Applied pressure [psi]';
        h1.Colormap = parula;
        h1.ColorbarVisible = 'off';
        h1.MissingDataColor = [0.5, 0.5, 0.5];
        h1.GridVisible = 'off';

        for s=1:length(sigma)
            subplot(1,1+length(sigma),1+s)
            h2=heatmap(cd,fliplr(delp),flipud(pO_values(:,:,s)));
            h2.Title = ['\fontsize{15} Osmotic Pressure at \sigma = ',num2str(sigma(s))];
            h2.XLabel = '\fontsize{15} Dialysate concentration [mM]';
            h2.YLabel = '\fontsize{15} Applied pressure [psi]';
            h2.Colormap = parula;
            h2.ColorbarVisible = 'off';
            h2.MissingDataColor = [0.5, 0.5, 0.5];
            h2.GridVisible = 'off';
        end

        if length(cd) <  11
            set(gcf,'Units','Inches','Position',[0.01 0.01 8*(length(sigma)+1) 5])
        else
            set(gcf,'Units','Inches','Position',[0.01 0.01 16*(length(sigma)+1) 10])
        end

        plotname = ['osmotic_sensitivity_heatmap-',model_stru.filenamestr,'.png'];
        saveas(gcf,plotname)
        movefile(plotname,'sigma_sensitivity')
    end

    figure(3)
    [X,Y] = meshgrid(cd,delp);
    yy = Y(:);
    colorstring = 'rbg';

    for s=1:length(sigma)
        pi_os = pO_values(:,:,s);
        pp = pi_os(:);
        subplot(length(sigma),3,1+3*(s-1))
        hold on
        mfig(s) = plot(yy./pp,range_mV(:),'o','Color',colorstring(s));
        xlabel('Difference in Mass Predictions','FontSize',15)
        ylabel('Applied pressure / Osmotic pressure','FontSize',15);

        subplot(length(sigma),3,2+3*(s-1))
        hold on
        cffig(s) = plot(yy./pp,range_cF(:),'o','Color',colorstring(s));
        xlabel('Difference in Retentate Predictions','FontSize',15)
        ylabel('Applied pressure / Osmotic pressure','FontSize',15);  

        subplot(length(sigma),3,3+3*(s-1))
        hold on
        cpfig(s) = plot(yy./pp,range_cH(:),'o','Color',colorstring(s));    
        xlabel('Difference in Permeate Predictions','FontSize',15)
        ylabel('Applied pressure / Osmotic pressure','FontSize',15);    

    end

    set(gcf,'Units','Inches','Position',[0.01 0.01 24 5*(length(sigma))])

    plotname = ['delP-diff-',model_stru.filenamestr,'.png'];
    saveas(gcf,plotname)
    movefile(plotname,'sigma_sensitivity')
end
end