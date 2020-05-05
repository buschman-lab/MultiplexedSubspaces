function FitMotifs(fn,save_fn,chunk,varargin)
%Camden MacDowell - timeless
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
end
rng('default'); %for reproducibility

%parse optional inputs
gp = general_params;
gp = ParseOptionalInputs(gp,varargin);

%load the training data
temp = load(fn,'data_train','data_test');
data_train = squeeze(temp.data_train(:,chunk,:));
data_test = squeeze(temp.data_test(:,chunk,:));

%fit lambda. don't need as many iterations 
lambda = FitLambda(data_train,gp);

%fit motifs 
fprintf('\n\n\nFitting Motifs')
[w, h, cost] = seqNMF(data_train, ...    
    'K', gp.K, 'L',gp.L, 'lambda',lambda,...        
    'showPlot', 0, 'maxiter',gp.maxiter,'tolerance',gp.tolerance,...
    'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
    'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',gp.shift);

%gather data
stats_train = CNMFStats(w,h,data_train);
stats_train.cost = cost;

%cross validate
fprintf('\n\n\nCrossvalidating Motifs')
[w_test, h_test] = seqNMF(data_test, ...    
    'K', gp.K, 'L',gp.L, 'lambda',lambda,...        
    'showPlot', 0, 'maxiter',gp.maxiter,'tolerance',gp.tolerance,...
    'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
    'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',gp.shift,...
    'W_init',w,'W_fixed',1);

%gather data
stats_test = CNMFStats(w_test,h_test,data_test);

%save off the data in the scratch directory and the nanpxs
fprintf('\n\tSaving data')

%save off the information in the scratch directory
save(save_fn,'stats_test','stats_train','w','-v7.3')

fprintf('\n\tDONE')

end
















