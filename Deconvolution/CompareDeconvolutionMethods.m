%Compare_DeconvolutionMethods
%Camden MacDowell - timeless

%grab recordings
% dff_list = GrabFiles(%
dff_list = {'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_08_2021_1dff_combined.mat'};
% spike_opts_list= GrabFiles(%
spike_opts_list = {'G:\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat'};

%load and parse data
[dff,st] = CompileData_deconvolution(dff_list,spike_opts_list);

stats = cell(1,4);
%narx neural net deconvolution
[ypred,y] = Deconvolve_NeuralNetwork(dff,st);
stats{1} = deconvolutionfitStats(ypred,y);

%glm kernel deconvolution
[ypred,y] = Deconvolve_GLM(dff,st,'direct');
stats{2} = deconvolutionfitStats(ypred,y);

%lucy richardson glm deconvolution
[ypred,y] = Deconvolve_GLM(dff,st,'lucyrichardson');
stats{3} = deconvolutionfitStats(ypred,y);

%lucy richardson GCAMP deconvolution
[ypred,y] = Deconvolve_GLM(dff,st,'lucyrichardson_gcamp');
stats{4} = deconvolutionfitStats(ypred,y);

%% make plots comparing the methods
%compare correlation to withheld data
figure; hold on; 
data = arrayfun(@(n) diag(stats{n}.rho), 1:numel(stats),'UniformOutput',0);
data = [data{:}];
bar(nanmean(data,1))

%compare generalizability across areas
figure; hold on; 
idx = eye(size(stats{1}.err));
data = arrayfun(@(n) stats{n}.err(~idx), 1:numel(stats),'UniformOutput',0);
data = [data{:}];
bar(nanmean(data,1))



