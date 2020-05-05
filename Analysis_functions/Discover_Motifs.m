function [w, cost, loadings, power, numFactors,...
    Xhat, Residuals, W, H, ExpVar_all, ExpVar_frame] = Discover_Motifs(X,varargin)

%THIS IS NOW DEPRECIATED AS OF 4/20/2020
%Set options
opts = general_params;
opts = ParseOptionalInputs(opts,varargin); 

%for reproducibility
rng('default'); 


%Run seqNMF for short snippet
w = seqNMF(X(:,1:10000), ...    
  'K', opts.K, 'L',opts.L, 'lambda',opts.lambda,...        
  'showPlot', 0, 'maxiter', 50,'tolerance',opts.tolerance,...
  'SortFactors',0,'lambdaL1H',opts.lambdaL1H,...
  'lambdaOrthoH',opts.lambdaOrthoH,'useWupdate',0,'Shift',opts.shift);
toc

tic
%Run seqNMF
[w, h, cost, loadings, power] = seqNMF(X, ...    
  'K', opts.K, 'L',opts.L, 'lambda',opts.lambda,...        
  'showPlot', 0, 'maxiter',2,'tolerance',opts.tolerance,...
  'SortFactors',0,'lambdaL1H',opts.lambdaL1H,...
  'lambdaOrthoH',opts.lambdaOrthoH,'useWupdate',0,'Shift',opts.shift);
toc
%Gather and organize seqNMF outputs 
[w, cost, loadings, power, numFactors, Xhat, Residuals, W, H] = ...
    CollectNMFOutputs(w, h, cost, loadings, power, nanpxs, opts.block, X);

%Get Explained variance
ExpVar_all = CalculateExplainedVariance(X,Residuals);
ExpVar_frame = CalculateExplainedVarianceFrameWise(X,Residuals);





end %end function loop








