function [stats_test, stats_train, w, h, gp, fh] = FitMotifs(data_train,data_test,varargin)
%Camden MacDowell - timeless

rng('default'); %for reproducibility

%parse optional inputs
gp = general_params;
gp = ParseOptionalInputs(gp,varargin);

%Optionally Fit Lambda (slow, but best practice).
if gp.lambda == -1 %fit lambda
    [lambda, fh] = FitLambda(data_train,gp,1);    
else %use input value 
    lambda = gp.lambda; 
end

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
stats_train.lambda = lambda;

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

%option to only keep motifs that contribute to cross validation and also to
%do the reverse fit


end
















