function GenerateExampleTraces()
%Camden MacDowell - timeless
%from CompareDeconvolutionMethods, but for generating data for figures
%provide params and training options to recreate analysis

rng('default')

%% Set parameters
norm_method = 'std';
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\';

%% Compile data
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\restingstate_processed_fn.mat','dff_list','spike_opts_list') 

%compile imaging data
[dff,st] = CompileData_deconvolution(dff_list,spike_opts_list,params);
close all;      

%normalize, split train/test and train
n = floor(size(dff{1},1)*3/4);
switch norm_method
    case 'mapminmax'
        dff_train = cellfun(@(x) mapminmax(x(1:n,:)',0,1)',dff,'UniformOutput',0);
        dff_test = cellfun(@(x) mapminmax(x(n+1:end,:)',0,1)',dff,'UniformOutput',0);            
        st_train = cellfun(@(x) mapminmax(x(1:n,:)',0,1)',st,'UniformOutput',0);
        st_test = cellfun(@(x) mapminmax(x(n+1:end,:)',0,1)',st,'UniformOutput',0);   
        %train all methods - optionally can adjust num neurons to mapmin
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',1000,params.bindata);                 
    case 'std'
        dff_train = cellfun(@(x) x(1:n,:)./std(x(1:n,:)),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:)./std(x(n+1:end,:)),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) x(1:n,:)./std(x(1:n,:)),st,'UniformOutput',0);
        st_test = cellfun(@(x) x(n+1:end,:)./std(x(n+1:end,:)),st,'UniformOutput',0);  
        %train all methods - optionally can adjust num neurons to std
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',100,params.bindata);                
    case 'none'
        dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
        dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
        st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
        st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
        %train all methods using recordings (10k sum fr across entire probe over entire rec) 
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
    otherwise
        error('unknown normalization method');
end  %switch end
% 
%loop through each recording and probe and deconvolution type
n_rec = numel(dff);
n_probe = size(dff{1},2);
ypred_ff = cell(1,n_rec);
ypred_glm = cell(1,n_rec);
ypred_lr = cell(1,n_rec);
ypred_none = cell(1,n_rec);
for cur_rec = 1:n_rec
    fprintf('\n\t working on rec %d or %d',cur_rec,n_rec);
    for cur_probe = 1:n_probe
       fprintf('\n working on probe %d or %d',cur_probe,n_probe)
       ypred_ff{cur_rec}(:,cur_probe) = DeconvolveData(dff_test{cur_rec}(:,cur_probe),'feedforward',trained_opts{cur_rec}(:,cur_probe));
       ypred_glm{cur_rec}(:,cur_probe) = DeconvolveData(dff_test{cur_rec}(:,cur_probe),'glm',trained_opts{cur_rec}(:,cur_probe));
       ypred_lr{cur_rec}(:,cur_probe) = DeconvolveData(dff_test{cur_rec}(:,cur_probe),'lr_gcamp',trained_opts{cur_rec}(:,cur_probe));
       ypred_none{cur_rec}(:,cur_probe) = DeconvolveData(dff_test{cur_rec}(:,cur_probe),'none',trained_opts{cur_rec}(:,cur_probe));
    end %probe
end %rec
save([savedir filesep sprintf('deconvolved_traces_withheld_data%s.mat',norm_method)],'params','st_test','ypred_ff','ypred_glm','ypred_lr','ypred_none');

%on the trained data
n_rec = numel(dff);
n_probe = size(dff{1},2);
ypred_ff = cell(1,n_rec);
ypred_glm = cell(1,n_rec);
ypred_lr = cell(1,n_rec);
ypred_none = cell(1,n_rec);
for cur_rec = 1:n_rec
    fprintf('\n\t working on rec %d or %d',cur_rec,n_rec);
    for cur_probe = 1:n_probe
       fprintf('\n working on probe %d or %d',cur_probe,n_probe)
       ypred_ff{cur_rec}(:,cur_probe) = DeconvolveData(dff_train{cur_rec}(:,cur_probe),'feedforward',trained_opts{cur_rec}(:,cur_probe));
       ypred_glm{cur_rec}(:,cur_probe) = DeconvolveData(dff_train{cur_rec}(:,cur_probe),'glm',trained_opts{cur_rec}(:,cur_probe));
       ypred_lr{cur_rec}(:,cur_probe) = DeconvolveData(dff_train{cur_rec}(:,cur_probe),'lr_gcamp',trained_opts{cur_rec}(:,cur_probe));
       ypred_none{cur_rec}(:,cur_probe) = DeconvolveData(dff_train{cur_rec}(:,cur_probe),'none',trained_opts{cur_rec}(:,cur_probe));
    end %probe
end %rec
st_test = st_train;
% save([savedir filesep sprintf('deconvolved_traces_trained_data%s.mat',norm_method)],'params','st_test','ypred_ff','ypred_glm','ypred_lr','ypred_none');
save([savedir filesep sprintf('deconvolved_traces_trained_data_10pvalidationAndGLMintercept%s.mat',norm_method)],'params','st_test','ypred_ff','ypred_glm');

end %function end


