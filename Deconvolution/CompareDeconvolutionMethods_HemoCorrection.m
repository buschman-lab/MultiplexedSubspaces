function CompareDeconvolutionMethods_HemoCorrection(norm_method)
%Camden MacDowell - timeless
%Compares the efficacy and generalizatbily of multiple deconovlution
%methods
%each comparison is assigned a separet 'type' to allow for easy spock

%addpaths for spock
if ispc
   dff_corrected = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse360_RestingState_10_17_2021_1dff_combined.mat';
   dff_uncorrected = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse360_RestingState_10_17_2021_1dff_combined_dff_uncorrected.mat';
   spike_opts_list = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse360_10_17_2021_RestingStateHemo_g0\ap_opts.mat';
else
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))  
   dff_corrected = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/Mouse360_RestingState_10_17_2021_1dff_combined.mat';
   dff_uncorrected = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/Mouse360_RestingState_10_17_2021_1dff_combined_dff_uncorrected.mat';
   spike_opts_list = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/ProcessedEphys/catgt_Mouse360_10_17_2021_RestingStateHemo_g0/ap_opts.mat';
end   

%compile imaging data
fprintf('\n\t Compiling Imaging Data')
params.bindata = 0; %temporally bin the imaging data. (don't use for hemo)
params.radius = 2; %pixel radius around probe tip    
params.offset = {[0,0]}; %optional moves center of DFF circle for probes at angles. [x,y]
[dff,~] = CompileData_deconvolution({dff_corrected,dff_uncorrected},[],params); %first is correct, second is uncorrected

% dff{2} = filterstack(dff{2}, 15, [0.1 4], 'lowpass', 1, 0);

%Compile ephys data (time bin to match corrected FPS). If depths have not
%been added then run AddVerticalDepthsToSpikeOpts. before continuing
fprintf('\n\t Compiling Ephys Data')
params.bindata = 1; %now time bin the ephys to match
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe
[~,st,n_neurons] = CompileData_deconvolution([],{spike_opts_list},params);
n_neurons = repmat(n_neurons,1,numel(dff));
st = repmat(st,1,numel(dff));
  
%save directory
if ispc
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\hemocorrected\';
else
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Deconvolution/hemocorrected/';
end
if ~exist(savedir,'dir'); mkdir(savedir); end

%normalize
switch norm_method
    case 'mapminmax'
        st = cellfun(@(x) mapminmax(x',0,1)',st,'UniformOutput',0);                
    case 'std'     
        st = cellfun(@(x) x./std(x),st,'UniformOutput',0);                
    case 'none'
        %do nothing
    case 'mean'
        for i = 1:numel(st)
           st{i} = st{i}./n_neurons(i,:);
        end                   
    otherwise
        error('unknown normalization method');
end

%%
%split into 6 ~10 minute chunks and train/test
t = 200:size(st{1},1); %burn in period (in frames) for hemocorrection to kick in
num_chunks = 10;
chunk_dur = floor(numel(t)/num_chunks);

t = arrayfun(@(x) t( ((x-1)*chunk_dur)+1:(x*chunk_dur)), 1:num_chunks,'UniformOutput',0);

stats_lr = cell(numel(t),numel(st));
stats_glm = cell(numel(t),numel(st));
stats_fNN = cell(numel(t),numel(st));
stats_none = cell(numel(t),numel(st));
stats_lr_train = cell(numel(t),numel(st));
stats_glm_train = cell(numel(t),numel(st));
stats_fNN_train = cell(numel(t),numel(st));
stats_none_train = cell(numel(t),numel(st));
grp = cell(numel(t),numel(st));
for i = 1:numel(t)
    %train on all except the given chunk 
    idx = 1:num_chunks; 
    idx = idx(~ismember(idx,i)); %remove the train setting
    dff_train = cellfun(@(x) x([t{idx}],1), dff,'UniformOutput',0);
    st_train = cellfun(@(x) x([t{idx}],1), st,'UniformOutput',0);    
    trained_opts = Deconvolve_Train(dff_train,st_train,'all',0,params.bindata);

%     %train on one chunk and fit to the next
%     dff_train = cellfun(@(x) x(t{i},1), dff,'UniformOutput',0);
%     st_train = cellfun(@(x) x(t{i},1), st,'UniformOutput',0);    
%     trained_opts = Deconvolve_Train(dff_train,st_train,'all',0,params.bindata);
%     
%     if i == numel(t)
%         dff_test = cellfun(@(x) x(t{1},1), dff,'UniformOutput',0);
%         st_test = cellfun(@(x) x(t{1},1), st,'UniformOutput',0);  
%     else
%         dff_test = cellfun(@(x) x(t{i+1},1), dff,'UniformOutput',0);
%         st_test = cellfun(@(x) x(t{i+1},1), st,'UniformOutput',0);  
%     end
    %fit to the withheld chunk
    dff_test = cellfun(@(x) x(t{i},1), dff,'UniformOutput',0);
    st_test = cellfun(@(x) x(t{i},1), st,'UniformOutput',0);    
    
    %evaluate on the trained sections
    lrpred={}; glmpred={}; fNNpred={}; nonepred={}; 
    for j = 1:numel(st)
        grp{i,j} = j;
        [stats_lr{i,j},lrpred{j}] = Deconvolve_Test(dff_test{j},st_test{j},'lr_gcamp',trained_opts{j});        
        [stats_glm{i,j},glmpred{j}] = Deconvolve_Test(dff_test{j},st_test{j},'glm',trained_opts{j},1);
        [stats_fNN{i,j},fNNpred{j}] = Deconvolve_Test(dff_test{j},st_test{j},'feedforward',trained_opts{j});
        [stats_none{i,j},nonepred{j}] = Deconvolve_Test(dff_test{j},st_test{j},'none',trained_opts{j});
    end
    sttrue = st_test{j};
    
    lrpred_train={}; glmpred_train={}; fNNpred_train={}; nonepred_train={}; 
    for j = 1:numel(st)
        grp{i,j} = j;
        [stats_lr_train{i,j},lrpred_train{j}] = Deconvolve_Test(dff_train{j},st_train{j},'lr_gcamp',trained_opts{j});        
        [stats_glm_train{i,j},glmpred_train{j}] = Deconvolve_Test(dff_train{j},st_train{j},'glm',trained_opts{j},1);
        [stats_fNN_train{i,j},fNNpred_train{j}] = Deconvolve_Test(dff_train{j},st_train{j},'feedforward',trained_opts{j});
        [stats_none_train{i,j},nonepred_train{j}] = Deconvolve_Test(dff_train{j},st_train{j},'none',trained_opts{j});
    end    
    sttrue_train = st_train{j};
end


grp = cat(1,grp{:});
stats_lr_train = cat(1,stats_lr_train{:});
stats_glm_train = cat(1,stats_glm_train{:});
stats_fNN_train = cat(1,stats_fNN_train{:});
stats_none_train = cat(1,stats_none_train{:});


fprintf('\n\tsaving')
%save off        
save([savedir filesep sprintf('hemocomparestats%s.mat',norm_method)],'grp','stats_lr','stats_glm','stats_fNN',...
    'trained_opts','params','stats_none','norm_method','n_neurons');   
save([savedir filesep sprintf('hemocompare_exampletraces%s.mat',norm_method)],'lrpred','glmpred','fNNpred','nonepred','sttrue');   

save([savedir filesep sprintf('hemocomparestats_train%s.mat',norm_method)],'grp','stats_lr_train','stats_glm_train','stats_fNN_train',...
    'trained_opts','params','stats_none_train','norm_method','n_neurons');   
save([savedir filesep sprintf('hemocompare_exampletraces_train%s.mat',norm_method)],'lrpred_train','glmpred_train','fNNpred_train','nonepred_train','sttrue_train'); 

end






