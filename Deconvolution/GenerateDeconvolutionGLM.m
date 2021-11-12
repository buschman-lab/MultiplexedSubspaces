function GenerateDeconvolutionGLM(norm_method,allrecs)
%Camden MacDowell - timeless
%See CompareDeconvolutionMethods for alternatives and paramter choices
%Fits a feedforward neural netowrk to the entire data set (all sites) for a
%recording using that recordings marked probe locations (preprocessing)
%Generates the FitNN used in ProcessAndSplitData. Takes a long time to run
%so reccomended to run independently prior to motif discovery
%@input type: if allrecs=1 (def) trains across all recs and sites in
%restinstats_processed_fn.mat. Otherwise trains within each recording

%addpaths for spock
if ispc
   load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list','spike_opts_list') 
else
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
   load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/restingstate_processed_fn.mat','dff_list_bucket','spike_opts_list_bucket')   
   dff_list = dff_list_bucket;
   spike_opts_list = spike_opts_list_bucket;
end   

%compile imaging data
fprintf('\n\t Compiling Imaging Data')
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
[dff,~] = CompileData_deconvolution(dff_list,[],params);
close all; 

%save directory where it will save off the neural network
if ispc
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\';
else
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/';    
end

params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 600]; %depth from surface of probe
%compile spiking data
[~,st,n_neurons] = CompileData_deconvolution([],spike_opts_list,params);   

%split train/test (additional train,validation,test occurs within train)
switch norm_method                
    case 'std'      
        st = cellfun(@(x) x(1:end,:)./std(x(1:end,:)),st,'UniformOutput',0);                
    case 'mean'
        for i = 1:numel(st)
           st{i} = st{i}./n_neurons(i,:);
        end              
    otherwise
        error('unknown normalization method');
end

%evaulate the fit
dff_train = dff; 
st_train = st;

if allrecs %train across everything       
    %combine training data xsites (Deconvolve_Train uses last 20% for
    %validation so make sure the last 20% for each rec is split out.
    n = size(dff_train{1},1);
    fprintf('\n\t n = %d',n);
    [trainInd,valInd] = divideblock(n,0.8, 0.2);
    temp_train_rec = cellfun(@(x) x(trainInd,:),dff_train,'UniformOutput',0); %get training and flatten across probes
    temp_train_rec = cellfun(@(x) x(:),temp_train_rec,'UniformOutput',0); 
    temp_train_rec = cat(1,temp_train_rec{:});
    temp_val_rec = cellfun(@(x) x(valInd,:),dff_train,'UniformOutput',0); %get val and flatten across probes
    temp_val_rec = cellfun(@(x) x(:),temp_val_rec,'UniformOutput',0); 
    temp_val_rec = cat(1,temp_val_rec{:});
    %combine across everyting 
    dff_train = {cat(1,temp_train_rec,temp_val_rec)};

    %same for the spiking data
    temp_train_rec = cellfun(@(x) x(trainInd,:),st_train,'UniformOutput',0); %get training and flatten across probes
    temp_train_rec = cellfun(@(x) x(:),temp_train_rec,'UniformOutput',0);
    temp_train_rec = cat(1,temp_train_rec{:});
    temp_val_rec = cellfun(@(x) x(valInd,:),st_train,'UniformOutput',0); %get val and flatten across probes
    temp_val_rec = cellfun(@(x) x(:),temp_val_rec,'UniformOutput',0); 
    temp_val_rec = cat(1,temp_val_rec{:});

    %combine across everyting 
    st_train = {cat(1,temp_train_rec,temp_val_rec)};      

    %train all methods
    trained_opts = Deconvolve_Train(dff_train,st_train,'glm',500,params.bindata);    
        
    %Train index is 1 rec and 1 probe but test index is all 6 rec and 4 probes
    probe_idx = repmat(1:4,numel(dff),1)';
    rec_idx = repmat(1:numel(dff),4,1);
    test_idx = [rec_idx(:),probe_idx(:)];
    train_idx = ones(size(test_idx));      

    %get list of not enough spiking
    badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

    %get the fit stats
    stats = cell(size(train_idx,1),1);   
    bad_probe_idx=[];
    for i = 1:size(train_idx,1)    
        fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
        %if either of the probes are crappy, skip
        if badprobe{train_idx(i,1)}(train_idx(i,2))==1
            %leave variables blank
            bad_probe_idx = [bad_probe_idx,i];
        else
            stats{i} = Deconvolve_Test(dff{test_idx(i,1)}(:,test_idx(i,2)),st{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));                         
        end
    end %rec    
        
else
    %combine training data across sites
    for i = 1:numel(dff_train)
        dff_train{i} = cat(1,dff_train{i}(:));
        st_train{i} = cat(1,st_train{i}(:));
    end

    %train all methods
    trained_opts = Deconvolve_Train(dff_train,st_train,'glm',500,params.bindata);    

    %get list of not enough spiking
    badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

    %Train index is rec and 1 (probe) but test index is rec and 4
    probe_idx = repmat(1:4,numel(dff_train),1)';
    rec_idx = repmat(1:numel(dff),4,1);
    test_idx = [rec_idx(:),probe_idx(:)];
    train_idx = [test_idx(:,1),ones(size(test_idx,1),1)];

    stats = cell(size(train_idx,1),1);   
    bad_probe_idx=[];
    for i = 1:size(train_idx,1)    
        fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
        %if either of the probes are crappy, skip
        if badprobe{train_idx(i,1)}(train_idx(i,2))==1
            %leave variables blank
            bad_probe_idx=[bad_probe_idx,i];
        else
            stats{i} = Deconvolve_Test(dff{test_idx(i,1)}(:,test_idx(i,2)),st{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
        end
    end %rec

end
fprintf('\n\tsaving')
%save off        
save([savedir filesep 'fit_GLM'],'stats','bad_probe_idx',...
    'badprobe','trained_opts','params','train_idx','test_idx','norm_method','n_neurons');          

end %function