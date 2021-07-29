function [dff_probe,st_probe] = CompileData_deconvolution(dff_list,spike_opts_list, params)
%Camden MacDowell - timeless
%takes lits of files and folders and pcompile the data needed for
%deconoltuoin

%adjustable parameters
if nargin <3
    params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
    params.depth = 1000; %depth from surface of probe
    params.radius = 2; %pixel radius around probe tip    
    params.offset = {[1,0],[0,0],[-1,1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
    params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
end

%loop through files
for cur_file = 1:numel(dff_list)
    fprintf('\n\t loading data file %d of %d',cur_file, numel(dff_list));
    spike_opts_fn = spike_opts_list{cur_file};
    dff_fn = dff_list{cur_file};
    
    %load the imaging data
    data = load(dff_fn);
    dff = data.dff; 
    [x,y,z] = size(dff); 
    [dff,nanpxs] = conditionDffMat(dff); 
    coords = data.probe_coords; 
    opts = data.opts; 
    clear data %free up space
    
    %load ephys opts
    spike_opts = load(spike_opts_fn);
    spike_opts = spike_opts.opts;         

    %parse the radius around the probes
    dff_probe = NaN(z,numel(coords));
    for cur_probe = 1:numel(coords)
        %create mask of radius around probe tip
        temp = coords{cur_probe};
        mask = zeros(x,y);  
        temp(1,1)=temp(1,1)+params.offset{cur_probe}(1);
        temp(1,2)=temp(1,2)+params.offset{cur_probe}(2);
        mask(temp(1,2)-params.radius:temp(1,2)+params.radius,temp(1,1)-params.radius:temp(1,1)+params.radius)=1;

        %flatten
        mask(nanpxs)=[]; 

        dff_probe(:,cur_probe) = nanmean(dff(:,mask==1),2);    
    end
    %free up space
    clear dff; 
    
    %parse spikes
    st_probe = NaN(z,numel(coords));
    for cur_probe = 1:numel(coords)
        %load the ephys data
        spikes = load([spike_opts.nidaq_path,sprintf('AP_Probe%d.mat',cur_probe)]); 
        %grab top x mm of spike IDs
        if params.mua == 1 %include labeleds muas
           idx = 1:numel(spikes.clust_info.depth);       
        else
           idx = sum(ismember(spikes.clust_info.group,'good '),2)==5;       
        end
        IDs = spikes.clust_info.id(spikes.clust_info.depth(idx)<params.depth);     
        st = spikes.clust_info.spike_times(ismember(spikes.clust_info.spike_cluster,IDs)); 

        %bin to the framerate of the imaging
        fileID = fopen([spike_opts.nidaq_path,'CameraFrameFrontEdgeTimes.txt'],'r');
        formatSpec = '%f';   
        im_times = fscanf(fileID,formatSpec);

        %mua per bin   
        st = histcounts(st,im_times);
        st = st/std(st); %std normalize

        %pad with single sample at end to match for ease
        st_probe(:,cur_probe) = [st,0]';
    end




end %function end




















