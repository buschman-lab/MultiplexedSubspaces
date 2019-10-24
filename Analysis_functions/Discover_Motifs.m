function [ExpVar_frame,ExpVar_all,Xhat,...
    Residuals,H,cost,loadings,power,numFactors,W,w,opts] = Discover_Motifs(varargin)

%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data','nanpxs');

fprintf('\nDiscovering Motifs Using Kval of %d \n',opts.K)

%Get current data block
X = data{opts.block};
fprintf('\tDiscovering Motifs for data group %d of %d...\n',opts.block,size(data,2))

%for reproducibility
rng('default'); 

%Run seqNMF
[w, h, cost, loadings, power] = seqNMF(X, ...    
  'K', opts.K, 'L',opts.L, 'lambda',opts.lambda,...        
  'showPlot', opts.showPlot, 'maxiter', opts.maxiter,'tolerance',opts.tolerance,...
  'SortFactors',1,'lambdaL1H',opts.lambdaL1H,...
  'lambdaOrthoH',opts.lambdaOrthoH,'useWupdate',1,'Shift',opts.shift);

%Gather and organize seqNMF outputs 
[w, cost, loadings, power, numFactors, Xhat, Residuals, W, H] = ...
    CollectNMFOutputs(w, h, cost, loadings, power, nanpxs, opts.block, X);

%Get Explained variance
ExpVar_all = CalculateExplainedVariance(X,Residuals);
ExpVar_frame = CalculateExplainedVarianceFrameWise(X,Residuals);

end %end function loop








