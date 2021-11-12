function CompareDeconvolution_NormalizationWindowLength(block,norm_method)
%Camden MacDowell - timeless
%Compares performance of deconvolution methods across different durations 
%dff normalziation windows
%block is which preprocessed data to do it on (i.e. 6 recs x # of windows)

%addpaths for spock
if ispc
   load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list','spike_opts_list') 
   data_path = GrabFiles('\w*seconds.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Normalization'});
else
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
   load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/restingstate_processed_fn.mat','dff_list_bucket','spike_opts_list_bucket')   
   dff_list = dff_list_bucket; %list of the data to match with spike opts
   spike_opts_list = spike_opts_list_bucket;   
   data_path = GrabFiles('\w*seconds.mat',0,{'/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Normalization'});
end   

%get the dff's to use and match to the spikes options
[~,fn_header] = cellfun(@(x) fileparts(x), dff_list,'UniformOutput',0);
header = cellfun(@(x) erase(x,'dff_combined'), fn_header,'UniformOutput',0);

%get the ephys rec that matches the current data
data_fn = data_path{block};
idx = cellfun(@(x) ~isempty(regexp(data_fn,x, 'once')),header,'UniformOutput',1);%find the header with the same recording
spike_opts_list = spike_opts_list(idx==1);

%compile imaging data
fprintf('\n\t Compiling Imaging Data')
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
[dff,~] = CompileData_deconvolution({data_fn},[],params);
close all; 

%train within site 
fprintf('\n\t Comparing within rec/mouse across sites')
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe
%compile spiking data
[~,st,n_neurons] = CompileData_deconvolution([],spike_opts_list,params);   

%depending on the question, dff are not required to be the whole recording, so adjust if necessary
n = size(dff{1},1);
st = cellfun(@(x) x(1:n,:), st,'UniformOutput',0);

%save directory
if ispc
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Normalization\WindowSweepResults\';
else
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Normalization/WindowSweepResults/';
end
if ~exist(savedir,'dir'); mkdir(savedir); end

%split train/test (additional train,validation,test occurs within train)
n = floor(size(dff{1},1)*3/4);
switch norm_method
    case 'mapminmax'
        dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) mapminmax(x(1:n,:)',0,1)',st,'UniformOutput',0);
        st_test = cellfun(@(x) mapminmax(x(n+1:end,:)',0,1)',st,'UniformOutput',0);   
        %train all methods - optionally can adjust num neurons to mapmin
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',1000,params.bindata);                 
    case 'std'
        dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) x(1:n,:)./std(x(1:n,:)),st,'UniformOutput',0);
        st_test = cellfun(@(x) x(n+1:end,:)./std(x(n+1:end,:)),st,'UniformOutput',0);  
        %train all methods - optionally can adjust num neurons to std
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',2500,params.bindata);                
    case 'none'
        dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
        st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
        %train all methods using recordings (10k sum fr across entire probe over entire rec) 
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
    case 'mean'
        for i = 1:numel(st)
           st{i} = st{i}./n_neurons(i,:);
        end
        dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
        st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
        %train all methods using recordings (10k sum fr across entire probe over entire rec) 
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',0,params.bindata);                
    otherwise
        error('unknown normalization method');
end

% Test within each probe per recroding
probe_idx = repmat(1:4,numel(dff_test),1)';
rec_idx = repmat(1:numel(dff),4,1);
train_idx = [rec_idx(:),probe_idx(:)];
%test on the same data (last quarter of the recording: see 'n' above)
test_idx = train_idx;             

%get list of not enough spiking
badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

%get number of units per rec
numunits = cellfun(@(x) NumUnitsPerDepthBin(x,params.depth), spike_opts_list,'UniformOutput',0);

num_methods = 4;
stats = cell(size(train_idx,1),num_methods);
stats_train = cell(size(train_idx,1),num_methods);
nunits_test = zeros(size(train_idx,1),num_methods);
nunits_train = zeros(size(train_idx,1),num_methods);
bad_probe_idx=[];
for i = 1:size(test_idx,1)                    
    fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
    %if either of the probes are crappy, skip
    if badprobe{train_idx(i,1)}(train_idx(i,2))==1 || badprobe{test_idx(i,1)}(test_idx(i,2))==1
        %leave variables blank
        bad_probe_idx=[bad_probe_idx,i];
    else
        stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
        stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        nunits_train(i,1:num_methods) = numunits{train_idx(i,1)}(train_idx(i,2));
        nunits_test(i,1:num_methods) = numunits{test_idx(i,1)}(test_idx(i,2));
        stats_train{i,1} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats_train{i,2} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
        stats_train{i,3} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        stats_train{i,4} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));                
    end
end %rec
%concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
deconv_stats = cat(1,stats{:});
deconv_stats_train = cat(1,stats_train{:});
type = {'lr_gcamp','glm','feedforward','none'};      
fprintf('\n\tsaving')

%parse the information about the data
[~,header] = fileparts(data_fn);
dur = str2double(regexp(header,'((?<=window).*(?=_))','match'));
mouse = str2double(regexp(header,'((?<=Mouse).*(?=_R))','match')); 
recdate = (regexp(header,'((?<=NP_).*(?=_1_win))','match'));

%save off
save([savedir filesep sprintf('comparelengths%s_block%d.mat',norm_method,block)],'deconv_stats','nunits_train','bad_probe_idx','type',...
    'nunits_test','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method','deconv_stats_train','n_neurons',...
    'dur','mouse','recdate','header');

end %function end

























