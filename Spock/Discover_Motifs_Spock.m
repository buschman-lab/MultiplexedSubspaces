function Discover_Motifs_Spock(varargin)

%Set options for fitting
opts = SetCNMFOptions(varargin);

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data','nanpxs');

fprintf('\nRunninng BuildRepitoire Using Kval of %d \n',opts.K)

%Make savedirectory
if exist(opts.save_dir)==0
    fprintf('Making save directory');
    mkdir(opts.save_dir);
end

%Get current data block
X = data{opts.block};
fprintf('\tFinding sequences for data group %d of %d...\n',block,size(data,2))

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
    CollectNMFOutputs(w, h, cost, loadings, power, opts.nanpxs, opts.block);

%Get Explained variance
[ExpVar_frame, ExpVar_all] = CalculateExplainedVariance(X,Xhat,Residuals);

%Make a movie if you want
if opts.movie_flag
   FormatFittingForMovie(X,W,H,opts);
end

%Save off parameters
if opts.save_data
    filename = [opts.save_dir sprintf('block_%d.mat',opts.block)];
    save(filename,'ExpVar_frame','ExpVar_all','Xhat','Residuals','H','cost','loadings','power','numFactors','W','w','opts','-v7.3');
end

end








