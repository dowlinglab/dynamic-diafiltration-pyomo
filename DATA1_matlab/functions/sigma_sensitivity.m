function sigsen_stru = sigma_sensitivity(data_stru,model_stru,sigma,theta,LOUD)
% Perform a set of sigma sensitivity simulations 
% Arguments:   
%     data_stru - all experimental data
%     model_stru - partial model configuration

colorstring = 'rbg';

fit_result = [model_stru.filenamestr,'\fit_stru','.mat'];
if isempty(theta) && exist(fit_result,'file')
    fit_stru = load(fit_result).fit_stru;
    theta = fit_stru.x;
else
    theta = model_stru.theta0;
end

for s = 1:length(sigma)   
    theta(model_stru.para_ind.sigma(1):model_stru.para_ind.sigma(2)) = sigma(s);
    sigsen_stru(s).theta = theta;
    sigsen_stru(s).sim_stru = sim_model(theta, data_stru, model_stru);
end

if LOUD
    figure
    for s = 1:length(sigma)   
        sim_stru = sigsen_stru(s).sim_stru;
        simname = ['sim_stru-dat', num2str(data_stru.dataset), ' C_Fin', num2str(model_stru.initialization.C_F0), 'sig', num2str(sigma(s)), '.mat'];
        save(simname, 'sim_stru')
        lg_name{s} = ['sigma = ', num2str(sigma(s))];
        for i = 1:length(sim_stru)            
            % plot mass prediction
            subplot(1,3,1)
            hold on
            mfig(s) = plot(sim_stru(i).time, sim_stru(i).mV,'Color',colorstring(s),'LineWidth',2);

            % plot retentate concentration prediction 
            subplot(1,3,2)
            hold on
            for i = 1:length(sim_stru)
                cffig(s) = plot(sim_stru(i).time, sim_stru(i).cF,'Color',colorstring(s),'LineWidth',2);
            end   

            % plot permeate concentration prediction 
            subplot(1,3,3)
            hold on
            for i = 1:length(sim_stru)
                if isfield(sim_stru,'cH')
                    cpfig(s) = plot(sim_stru(i).time(:), sim_stru(i).cH,'Color',colorstring(s),'LineWidth',2);
                else
                    cpfig(s) = plot(sim_stru(i).time(2:end), sim_stru(i).cV(2:end),'Color',colorstring(s),'LineWidth',2);

                end
            end
        end
    end

    subplot(1,3,1)
    xlabel('Time [s]','FontSize',15)
    ylabel('Collected Permeate Mass [g]','FontSize',15);
    title({'\fontsize{15} Mass Predictions',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');
    set(gca,'FontSize',12)
    hold off

    subplot(1,3,2) 
    legend(cffig,lg_name,'FontSize',12,'Location','Best');
    xlabel('Time [s]','FontSize',15)
    ylabel(['Retentate Concentration of ',data_stru.data_config.namec,' [mmol/L]'],'FontSize',15);
    title({'\fontsize{15} Concentration Predictions',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');    
    set(gca,'FontSize',12)
    hold off

    subplot(1,3,3)
    xlabel('Time [s]','FontSize',15)
    ylabel(['Permeate Concentration of ',data_stru.data_config.namec,' [mmol/L]'],'FontSize',15);
    title({'\fontsize{15} Concentration Predictions',...
        '\fontsize{11}' model_stru.titlestr},'interpreter','tex');    
    set(gca,'FontSize',12)
    hold off

    set(gcf,'Units','Inches','Position',[0.01 0.01 24 5])
    figname = ['sigma_sensitivity-CF0_',num2str(model_stru.initialization.C_F0),'-',model_stru.filenamestr,'.png'];
    saveas(gcf,figname)
end
end