function FitTrainingData_junk(file_list,varargin)
%Camden MacDowell - timeless
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
end
rng('default'); %for reproducibility

%parse optional inputs
gp = general_params;
gp = ParseOptionalInputs(gp,varargin);

%load the training data
temp = load(file_list,'data_train','data_test');
data_train = temp.data_train;
data_test = temp.data_test;

gp.K = 25;
gp.L = 15;

fprintf('\nfitting motifs')
%find motifs 
[w, h, cost_list, loadings, power] = seqNMF(data_train, ...    
  'K', gp.K, 'L',gp.L, 'lambda',0.005,...        
  'showPlot', 0, 'maxiter',50,'tolerance',gp.tolerance,...
  'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
  'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',1);

seqNMF(data_train, ...    
  'K', gp.K, 'L',gp.L, 'lambda',gp.lambda,...        
  'showPlot', 1, 'maxiter',1,'tolerance',gp.tolerance,...
  'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
  'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',gp.shift,...
  'W_init',w,'H_init',h);

Xhat = helper.reconstruct(w,h);
Residuals = data_train-Xhat; 

%Get Explained variance
ExpVar_train = CalculateExplainedVariance(data_train,Residuals);


%Collect Information
p = (sum(data_train(:).^2)-sum((data_train(:)-Xhat(:)).^2))/sum(data_train(:).^2);  % fraction power explained by whole reconstruction

%Calc the correlation %remember camden, the rho is just the sqrt(variance)
rho = corr(data_train(:),Xhat(:));

[cost,reg,~] = helper.get_seqNMF_cost(data_train,w,h);
lambda = gp.lambda;

fprintf('\n cross-validating')
%cross validate on testing data
[w1, h1, ~, loadings_test, ~] = seqNMF(data_test, ...    
  'K', gp.K, 'L',gp.L, 'lambda',gp.lambda,...        
  'showPlot', gp.showPlot, 'maxiter',100,'tolerance',gp.tolerance,...
  'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
  'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',0,'Shift',gp.shift,...
  'W_init',w,'W_fixed',1);

Xhat = helper.reconstruct(w1,h1);
Residuals = data_test-Xhat;
ExpVar_test = CalculateExplainedVariance(data_train,Residuals);

%save off the 
fprintf('\n\tSaving data')
[~, fn_temp] = fileparts(file_list);  
fn = [gp.local_bucket gp.processing_intermediates fn_temp sprintf('_fit_lambda%d.mat',gp.lambda)];
if ~ispc
    fn = ConvertToBucketPath(fn);
end
save(fn,'w','h','cost','loadings','power','loadings_test','reg','cost','ExpVar_train','ExpVar_test','lambda','cost_list');
fprintf('\n\tDONE')



%% Plot the fits at different durations
COUNT = 1;
cost_list = {};
expvar_train = {};
expvar_test = {};
cost_list_norm = {};
label = {};
dur_val = [];
for dur = [15*[60,120,240,480,900]] %, 25900
    temp = load(sprintf('495_1_28_2020_1dff_combined_processed_fit_dur%d.mat',dur),'ExpVar_train','ExpVar_test','cost_list');
    cost_list{COUNT} = temp.cost_list;
    cost_list_norm{COUNT} = (temp.cost_list)/max(temp.cost_list);
    expvar_train{COUNT} = temp.ExpVar_train;
    expvar_test{COUNT} = temp.ExpVar_test;
    label{COUNT} = sprintf('Dur %d',dur);
    dur_val(COUNT) = dur;
    COUNT= COUNT+1;
end    

figure; hold on; 
title('RMSE per iteration across different durations')
cellfun(@plot, cost_list)
legend(label)
set(gca,'yscale','log','xlim',[0 300])

figure; hold on; 
title('RMSE per iteration across different durations (normalized)')
cellfun(@plot, cost_list_norm)
legend(label)
set(gca,'yscale','log','xlim',[0 300])

figure; hold on; 
title('explained variance of fit')
plot(dur_val,[expvar_train{:}],'linewidth',2,'Marker','o');
plot(dur_val,[expvar_test{:}],'linewidth',2,'Marker','o');
legend({'training','testing'});
xlabel('fit duration (frames)');
ylabel('pev');

handles = get(groot, 'Children');
saveCurFigs(handles,'-dpng',sprintf('troubleshootingfits',fn),pwd,0); %close all;












