function FitTrainingData(file_list,varargin)


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
data_train = temp.data_train(:,1:gp.dur);
data_test = temp.data_test(:,1:gp.dur);

fprintf('\nfitting motifs')
%find motifs 
[w, h, cost_list, loadings, power] = seqNMF(data_train, ...    
  'K', gp.K, 'L',gp.L, 'lambda',gp.lambda,...        
  'showPlot', 0, 'maxiter',gp.maxiter,'tolerance',gp.tolerance,...
  'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
  'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',gp.shift);

Xhat = helper.reconstruct(w,h);
Residuals = data_train-Xhat; 

%Get Explained variance
ExpVar_train = CalculateExplainedVariance(data_train,Residuals);

[cost,reg,~] = helper.get_seqNMF_cost(data_train,w,h);
lambda = gp.lambda;

fprintf('\n cross-validating')
%cross validate on testing data
[w1, h1, ~, loadings_test, ~] = seqNMF(data_test, ...    
  'K', gp.K, 'L',gp.L, 'lambda',gp.lambda,...        
  'showPlot', gp.showPlot, 'maxiter',gp.maxiter,'tolerance',gp.tolerance,...
  'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
  'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',0,'Shift',gp.shift,...
  'W_init',w,'W_fixed',1);

Xhat = helper.reconstruct(w1,h1);
Residuals = data_test-Xhat;
ExpVar_test = CalculateExplainedVariance(data_train,Residuals);

%save off the 
fprintf('\n\tSaving data')
[~, fn_temp] = fileparts(file_list);  
% fn = [gp.local_bucket gp.processing_intermediates fn_temp sprintf('_fit_lambda%d.mat',gp.lambda)];
fn = [gp.local_bucket gp.processing_intermediates fn_temp sprintf('_fit_dur%d.mat',gp.lambda)];
if ~ispc
    fn = ConvertToBucketPath(fn);
end
save(fn,'w','h','cost','loadings','power','loadings_test','reg','cost','ExpVar_train','ExpVar_test','lambda','cost_list');
fprintf('\n\tDONE')















