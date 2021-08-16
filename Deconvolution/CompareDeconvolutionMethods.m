function CompareDeconvolutionMethods(type,block)
%Camden MacDowell - timeless
%Compares the efficacy and generalizatbily of multiple deconovlution
%methods
%each comparison is assigned a separet 'type' to allow for easy spock

%addpaths for spock
if ispc
   load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\restingstate_processed_fn.mat','dff_list','spike_opts_list') 
else
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
   addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
   load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Deconvolution/restingstate_processed_fn.mat','dff_list_bucket','spike_opts_list_bucket')   
   dff_list = dff_list_bucket;
   spike_opts_list = spike_opts_list_bucket;
end   

%normalization method (so consistent across all) options =
%mapminmax, std, or none
norm_method = 'std';

%compile imaging data
fprintf('\n\t Compiling Imaging Data')
params.bindata = 0; %temporally bin the data?
params.radius = 2; %pixel radius around probe tip    
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
[dff,~] = CompileData_deconvolution(dff_list,[],params);
close all; 

%save directory
if ispc
    if params.bindata==1
        savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\timebinned\';
    else
        savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\';
    end
else
    if params.bindata==1
        savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Deconvolution/timebinned/';
    else
        savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Deconvolution/';
    end
end

    
switch type %different comparisons to run 
    case 1 %generalizability within recording, within site across depths
        fprintf('\n\t Comparing within rec/mouse/site across depths')
        depths = [100:200:900; 200:200:1000]';
        depths = cat(1,depths,[0 600],[600 1000]);
        
        %run for a given block (depth). Takes about 4-12 hours per depth so do this way
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = depths(block,:); %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);   

        %split train/test (additional train,validation,test occurs within train)
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
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',2500,params.bindata);                
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
                %train all methods using recordings (10k sum fr across entire probe over entire rec) 
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
            otherwise
                error('unknown normalization method');
        end          
        
        % Test within each probe per recroding
        probe_idx = repmat(1:4,numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),4,1);
        train_idx = [rec_idx(:),probe_idx(:)];
        %test on the same data (second have of the recording)
        test_idx = train_idx;             

        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        %get number of units per rec
        numunits = cellfun(@(x) NumUnitsPerDepthBin(x,params.depth), spike_opts_list,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);
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
                stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,5} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                nunits_train(i,1:num_methods) = numunits{train_idx(i,1)}(train_idx(i,2));
                nunits_test(i,1:num_methods) = numunits{test_idx(i,1)}(test_idx(i,2));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','lr_glm','glm','feedforward','none'};      
        fprintf('\n\tsaving')
        %save off
        save([savedir filesep sprintf('within_comparedepths%s%d.mat',norm_method,block)],'type','deconv_stats','nunits_train','bad_probe_idx',...
            'nunits_test','badprobe','trained_opts','params','depths','block','train_idx','test_idx','norm_method');
                
    case 2 %generalizability within recording, across sites same depth       
        fprintf('\n\t Comparing within rec/mouse across sites')
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = [0 600]; %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);          

        %split train/test (additional train,validation,test occurs within train)
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
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',2500,params.bindata);                
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
                %train all methods using recordings (10k sum fr across entire probe over entire rec) 
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
            otherwise
                error('unknown normalization method');
        end                 
        
        %each train will be used 3 times (other locations)
        temp = repmat(1:4,3,1);
        probe_idx = repmat(temp(:)',numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),numel(temp),1);
        train_idx = [rec_idx(:),probe_idx(:)];

        temp = repmat(1:4,4,1)';
        temp(logical(eye(size(temp)))) = []; %remove auto comparisons
        probe_idx = repmat(temp(:)',numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),numel(temp),1);
        test_idx = [rec_idx(:),probe_idx(:)];
        
        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);    
        bad_probe_idx=[];
        for i = 1:size(train_idx,1)    
            fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
            %if either of the probes are crappy, skip
            if badprobe{train_idx(i,1)}(train_idx(i,2))==1 || badprobe{test_idx(i,1)}(test_idx(i,2))==1
                bad_probe_idx=[bad_probe_idx,i];
                %leave variables blank
            else
                stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,5} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','lr_glm','glm','feedforward','none'};      

        fprintf('\n\tsaving')
        %save off        
        save([savedir filesep sprintf('within_xsites%s.mat',norm_method)],'type','deconv_stats','bad_probe_idx','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method');
        
    case 3 %generalizability across recording, same site/animal/depth
        fprintf('\n\t Comparing across rec, within mouse/sites')
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = [0 600]; %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);          

        %split train/test (additional train,validation,test occurs within train)
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
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',2500,params.bindata);                
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
                %train all methods using recordings (10k sum fr across entire probe over entire rec) 
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
            otherwise
                error('unknown normalization method');
        end
        %get recording ids 
        mousenum = MouseNumFromPath(dff_list,'Mouse_');
        unique_mice = unique(mousenum);

        %each train will be used onces
        probe_idx = repmat(1:4,numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),4,1);
        train_idx = [rec_idx(:),probe_idx(:)];

        %flip test days so diff than train days per mouse
        test_idx = train_idx; 
        for i = 1:numel(unique_mice)
            idx = find(mousenum==unique_mice(i));
            %flip the rec numbers
            for j = 1:numel(idx)
                test_idx(train_idx(:,1)==idx(1,j))=repmat(setdiff(idx,idx(j)),sum(train_idx(:,1)==idx(1,j)),1);
            end
        end
        
        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);   
        bad_probe_idx=[];
        for i = 1:size(train_idx,1)    
            fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
            %if either of the probes are crappy, skip
            if badprobe{train_idx(i,1)}(train_idx(i,2))==1 || badprobe{test_idx(i,1)}(test_idx(i,2))==1
                %leave variables blank
                bad_probe_idx=[bad_probe_idx,i];
            else
                stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,5} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','lr_glm','glm','feedforward','none'};      

        fprintf('\n\tsaving')
        %save off        
        save([savedir filesep sprintf('xrec_samesites%s.mat',norm_method)],'type','deconv_stats','bad_probe_idx','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method');
                
    case 4 %generalizability across animals, same site/depth
        fprintf('\n\t Comparing across animal, within sites')
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = [0 600]; %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);          

        %split train/test (additional train,validation,test occurs within train)
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
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',2500,params.bindata);                
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
                %train all methods using recordings (10k sum fr across entire probe over entire rec) 
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
            otherwise
                error('unknown normalization method');
        end 

        %get recording ids and make sure that recordings are in correct order
        mousenum = MouseNumFromPath(dff_list,'Mouse_');       
        assert(sum(mousenum([1,3,5])-mousenum([2,4,6]))==0,'mouse numbers not ordered correctly');

        % Same probe for each animals and then repeat, but with shifting
        % recording number
        probe_idx = repmat(1:4,numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),4,1);
        train_idx = [rec_idx(:),probe_idx(:)];
        test_idx = [];
        %now repeat such that anima1 1 d1 is tested on animal 2 day 1 and
        %animal 3 day 1
        for i = 1:3
            temp_test_idx = train_idx;
            if i == 1 %shift x2 to skip within animal
                temp_test_idx(:,1) = circshift(train_idx(:,1),8);
            else
                temp_test_idx(:,1) = circshift(train_idx(:,1),8+4*(i-1)); 
            end
            test_idx = cat(1,test_idx,temp_test_idx);
        end        
        %final round, at the single shift without self replicates
        train_idx=repmat(train_idx,3,1);
        pairings = [1,6;6,1;2,3;3,2;4,5;5,4];
        for i = 1:size(pairings,1)
            train_idx = cat(1,train_idx,[ones(4,1)*pairings(i,1),[1;2;3;4]]);
            test_idx = cat(1,test_idx,[ones(4,1)*pairings(i,2),[1;2;3;4]]);            
        end
          
        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);   
        bad_probe_idx=[];
        for i = 1:size(train_idx,1)    
            fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
            %if either of the probes are crappy, skip
            if badprobe{train_idx(i,1)}(train_idx(i,2))==1 || badprobe{test_idx(i,1)}(test_idx(i,2))==1
                %leave variables blank
                bad_probe_idx=[bad_probe_idx,i];
            else
                stats{i,1} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,3} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,4} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,5} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','lr_glm','glm','feedforward','none'};      

        fprintf('\n\tsaving')
        %save off        
        save([savedir filesep sprintf('xanimals_withinsites%s.mat',norm_method)],'type','deconv_stats','bad_probe_idx','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method');        
    
    case 5 %training data (add glm intercept and switch validation to 0.1 and test to 0 before running
        fprintf('\n\t Comparing within animal, on the training data')
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = [0 600]; %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);          

        %split train/test (additional train,validation,test occurs within train)
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
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',500,params.bindata);                
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
                %train all methods using recordings (10k sum fr across entire probe over entire rec) 
                trained_opts = Deconvolve_Train(dff_train,st_train,'all',10000,params.bindata);
            otherwise
                error('unknown normalization method');
        end 

        
        % Test within each probe per recroding
        probe_idx = repmat(1:4,numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),4,1);
        train_idx = [rec_idx(:),probe_idx(:)];
        %test on the same data (second have of the recording)
        test_idx = train_idx;   
          
        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);   
        bad_probe_idx=[];
        for i = 1:size(train_idx,1)    
            fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
            %if either of the probes are crappy, skip
            if badprobe{train_idx(i,1)}(train_idx(i,2))==1 || badprobe{test_idx(i,1)}(test_idx(i,2))==1
                %leave variables blank
                bad_probe_idx=[bad_probe_idx,i];
            else
                stats{i,1} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,2} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'lr_glm',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,3} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,4} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,5} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','lr_glm','glm','feedforward','none'};      

        fprintf('\n\tsaving')
        %save off        
        save([savedir filesep sprintf('training%s.mat',norm_method)],'type','deconv_stats','bad_probe_idx','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method');   
        
    case 6 %train across sites within recording
        fprintf('\n\t training within recording. Compare on both train and test data')
        params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
        params.depth = [0 600]; %depth from surface of probe
        %compile spiking data
        [~,st] = CompileData_deconvolution([],spike_opts_list,params);          

        %split train/test (additional train,validation,test occurs within train)
        n = floor(size(dff{1},1)*3/4);
        switch norm_method
            case 'mapminmax'
                dff_train = cellfun(@(x) mapminmax(x(1:n,:)',0,1)',dff,'UniformOutput',0);
                dff_test = cellfun(@(x) mapminmax(x(n+1:end,:)',0,1)',dff,'UniformOutput',0);            
                st_train = cellfun(@(x) mapminmax(x(1:n,:)',0,1)',st,'UniformOutput',0);
                st_test = cellfun(@(x) mapminmax(x(n+1:end,:)',0,1)',st,'UniformOutput',0);                  
            case 'std'
                dff_train = cellfun(@(x) x(1:n,:)./std(x(1:n,:)),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:)./std(x(n+1:end,:)),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:)./std(x(1:n,:)),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:)./std(x(n+1:end,:)),st,'UniformOutput',0);              
            case 'none'
                dff_train = cellfun(@(x) x(1:n,:),dff,'UniformOutput',0);
                dff_test = cellfun(@(x) x(n+1:end,:),dff,'UniformOutput',0);            
                st_train = cellfun(@(x) x(1:n,:),st,'UniformOutput',0);
                st_test = cellfun(@(x) x(n+1:end,:),st,'UniformOutput',0); 
            otherwise
                error('unknown normalization method');
        end 
        
        %combine training data xsites
        [trainInd,valInd,testInd] = divideblock(n,0.7, 0.2, 0.1);
        temp_train_rec = cellfun(@(x) x(trainInd,:),dff_train,'UniformOutput',0); %get training and flatten across probes
        temp_train_rec = cellfun(@(x) x(:),temp_train_rec,'UniformOutput',0); 
        temp_val_rec = cellfun(@(x) x(valInd,:),dff_train,'UniformOutput',0); %get val and flatten across probes
        temp_val_rec = cellfun(@(x) x(:),temp_val_rec,'UniformOutput',0); 
        temp_test_rec = cellfun(@(x) x(testInd,:),dff_train,'UniformOutput',0); %get test and flatten across probes
        temp_test_rec = cellfun(@(x) x(:),temp_test_rec,'UniformOutput',0); 
        for i = 1:numel(dff_train)
           dff_train{i} = cat(1,temp_train_rec{i},temp_val_rec{i},temp_test_rec{i});
        end
        
        %same for the spiking data
        temp_train_rec = cellfun(@(x) x(trainInd,:),st_train,'UniformOutput',0); %get training and flatten across probes
        temp_train_rec = cellfun(@(x) x(:),temp_train_rec,'UniformOutput',0); 
        temp_val_rec = cellfun(@(x) x(valInd,:),st_train,'UniformOutput',0); %get val and flatten across probes
        temp_val_rec = cellfun(@(x) x(:),temp_val_rec,'UniformOutput',0); 
        temp_test_rec = cellfun(@(x) x(testInd,:),st_train,'UniformOutput',0); %get test and flatten across probes
        temp_test_rec = cellfun(@(x) x(:),temp_test_rec,'UniformOutput',0); 
        for i = 1:numel(st_train)
           st_train{i} = cat(1,temp_train_rec{i},temp_val_rec{i},temp_test_rec{i});
        end          

        %train all methods
        trained_opts = Deconvolve_Train(dff_train,st_train,'all',500,params.bindata);    
        
        %Train index is rec and 1 (probe) but test index is rec and 4
        probe_idx = repmat(1:4,numel(dff_test),1)';
        rec_idx = repmat(1:numel(dff),4,1);
        test_idx = [rec_idx(:),probe_idx(:)];
        train_idx = [test_idx(:,1),ones(size(test_idx,1),1)];

        %get list of not enough spiking
        badprobe = cellfun(@(x) [x.insufactivity], trained_opts,'UniformOutput',0);

        num_methods = 6;
        stats = cell(size(train_idx,1),num_methods);   
        bad_probe_idx=[];
        for i = 1:size(train_idx,1)    
            fprintf('\t\n working on comparison %d of %d',i,size(train_idx,1));
            %if either of the probes are crappy, skip
            if badprobe{train_idx(i,1)}(train_idx(i,2))==1
                %leave variables blank
                bad_probe_idx=[bad_probe_idx,i];
            else
                stats{i,1} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'lr_gcamp',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,2} = Deconvolve_Test(dff_test{test_idx(i,1)}(:,test_idx(i,2)),st_test{test_idx(i,1)}(:,test_idx(i,2)),'glm',trained_opts{train_idx(i,1)}(train_idx(i,2))); 
                stats{i,3} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'feedforward',trained_opts{train_idx(i,1)}(train_idx(i,2)));
                stats{i,4} = Deconvolve_Test(dff_train{test_idx(i,1)}(:,test_idx(i,2)),st_train{test_idx(i,1)}(:,test_idx(i,2)),'none',trained_opts{train_idx(i,1)}(train_idx(i,2)));
            end
        end %rec
        %concatenate (results is all narx, then all lr_gcamp, then lr_glm, etc.
        deconv_stats = cat(1,stats{:});
        type = {'lr_gcamp','glm','feedforward','none'};      

        fprintf('\n\tsaving')
        %save off        
        save([savedir filesep sprintf('training%s.mat',norm_method)],'type','deconv_stats','bad_probe_idx','badprobe','trained_opts','params','block','train_idx','test_idx','norm_method');    
        
    case 7 %train across sites and animals and compare to withheld
            
        
    otherwise 
        error('unknown comparison')
end






