function Spock_CompareResolution(block,c_flag)
%Set path for spock to all repository (one step above cur dir)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end
if c_flag ==1 
    name = '431-10-17-2019_1dff_combined_uncorrected.mat';
else
    name = '431-10-17-2019_1dff_combined.mat';
end

%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,{'data_file_name',sprintf('/ProcessedData_Hemo/%s',name)}); 

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
stack = load([opts.bucket opts.base opts.data_file_name],'stack');
%for reproducibility
rng('default'); 

%get a random 2 minutes
for i = 1:block; rand(1); end %so random selection is depending on block number
start = randi(1,40000,1); 
dff = stack.stack(:,:,start:start+13*119);
clear stack; 

sb = {1,2,2};
sm = {0,0,1};
stat = struct();
for i = 1:3
    fprintf('working on round %d',i);
    if i ==1 
        [dff_norm_lin_array, nanpxs_array] = linearizeDff({dff},'sm_kern',[],'spatialbin',1,'filtband',[0.1 4],'nSTD',2);
    elseif i==2
        [dff_norm_lin_array, nanpxs_array] = linearizeDff({dff},'sm_kern',[],'spatialbin',2,'filtband',[0.1 4],'nSTD',2);
    elseif i==3
        [dff_norm_lin_array, nanpxs_array] = linearizeDff({dff},'sm_kern',[1 1 0.1],'spatialbin',2,'filtband',[0.1 4],'nSTD',2);
    end
    %linearize the data
    X = dff_norm_lin_array{1};
    nanpxs = nanpxs_array{1};

    fprintf('\nDiscovering Motifs Using Kval of %d \n',opts.K)

    %for reproducibility
    rng('default'); 

    %Run seqNMF
    [w, h, cost, loadings, power] = seqNMF(X, ...    
      'K', opts.K, 'L',opts.L, 'lambda',opts.lambda,...        
      'showPlot', 0, 'maxiter', opts.maxiter,'tolerance',opts.tolerance,...
      'SortFactors',1,'lambdaL1H',opts.lambdaL1H,...
      'lambdaOrthoH',opts.lambdaOrthoH,'useWupdate',1,'Shift',opts.shift);

    %Gather in case used GPU
    cost = {gather(cost)};
    loadings = {gather(loadings)};
    power = {gather(power)};
    w = gather(w);

    %Remove empty w or w with barely anything
    indempty = sum(sum(w>0,1),3)==0; % W is literally empty
    Wflat = sum(w,3); 
    indempty = indempty | (max(Wflat,[],1).^2> .5*sum(Wflat.^2,1)); % or one pixel has >50% of the power
    w(:,indempty,:) = []; % Delete factors that meet the above critera
    H = gather(h);
    H(indempty,:) = [];

    %Number of factors
    numFactors = {size(w,2)};

    %Get Residuals: 
    Xhat = helper.reconstruct(w,H);
    Residuals = X-Xhat; 

    %Get Explained variance
    stat(i).ExpVar_all = CalculateExplainedVariance(X,Residuals);
    stat(i).ExpVar_frame = CalculateExplainedVarianceFrameWise(X,Residuals);
    stat(i).numFactors = numFactors; 
    stat(i).spatialbin = sb{i};
    stat(i).smoothing = sm{i};

end

fprintf('saving...')
fn = [opts.bucket opts.base sprintf('/AnalyzedData_MesomappingManuscript_5_2019/CompareResolution/%s_block%d.mat',name,block)] ;
save(fn,'stat','-v7.3')
end %end function 







