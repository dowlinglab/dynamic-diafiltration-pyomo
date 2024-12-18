% Select dataset and modified initial feed concentration
% Single component experiments
% 501.1 filtration w/ conductivity measurements
% 511.12 diafiltration w/ conductivity measurements
Dat = 511.12;
% Sigma values evaluated at
sigma = [0.1,0.5,0.9];
% Parameters evaluated at
theta = [3.8994;0.29489;1];% NF90.5 diafiltration
% theta = [];

% Run simulation
data_stru = load(['data/data_stru-dataset', num2str(Dat), '.mat']).data_stru;
legacy_model_name = '201cvmv';

model_stru = model_config(legacy_model_name,data_stru);
model_stru.concpolar = true;
model_stru = update_model_config(data_stru,model_stru);
cF = readmatrix([data_stru.filename,'_cF.csv']);
if ~any(Dat==[501.12,511.12])
    model_stru.initialization.C_F0 = cF(1,2:end);
elseif Dat==511.12
    model_stru.initialization.C_F0 = cF(cF(:,1)==380,2:end);
end
sigsen_stru = sigma_sensitivity(data_stru, model_stru,sigma,theta,true);

plotname = ['sigma_sensitivity-CF0_',num2str(model_stru.initialization.C_F0),'-',model_stru.filenamestr,'.png'];
if ~isfolder('sigma_sensitivity')
    mkdir('sigma_sensitivity');
end
movefile(plotname,'sigma_sensitivity')









