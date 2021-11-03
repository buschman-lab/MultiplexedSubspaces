function [dff_probe,st_probe,n_neurons] = CompileData_deconvolution(dff_list,spike_opts_list,params,verbose)
%Camden MacDowell - timeless
%takes lits of files and folders and pcompile the data needed for
%deconvolution
%lead dff_list or spike_opts_list blank just compile one data type

%adjustable parameters
if nargin <3
    params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
    params.depth = [0 1000]; %depth from surface of probe
    params.radius = 2; %pixel radius around probe tip    
    params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
    params.bindata = 0; %1=temporally bin data (i.e. to frame rate 15 from 30)
end

if nargin <4
    verbose = 0;
end

%compile imaging data 
if ~isempty(dff_list)
    dff_probe = cell(1,numel(dff_list));
    for cur_file = 1:numel(dff_list)    
        fprintf('\n\t working on imaging data for file %d of %d',cur_file, numel(dff_list));
        dff_fn = dff_list{cur_file};
       
        %load imaging data
        data = load(dff_fn);
        dff = data.dff;
        if params.bindata ==1  
           dff = dff(1:2:end,:)+dff(2:2:end,:);             
           dff = conditionDffMat(dff,data.nanpxs,[]);
           dff = DenoisePCA(dff);
           [dff,nanpxs] = conditionDffMat(dff);        
        else
            nanpxs = data.nanpxs;
        end        
        coords = data.probe_coords;
        opts = data.opts;
        clear data;

        %get dimesions of the full data
        x = sqrt(size(dff,2)+numel(nanpxs));   
        z = size(dff,1);

        %parse the radius around the probes
        dff_probe{cur_file} = NaN(z,numel(coords));
        all_masks = cell(1,numel(coords));
        for cur_probe = 1:numel(coords)
            %create mask of radius around probe tip
            temp = coords{cur_probe};
            mask = zeros(x,x);  
            temp(1,1)=temp(1,1)+params.offset{cur_probe}(1);
            temp(1,2)=temp(1,2)+params.offset{cur_probe}(2);
            mask(temp(1,2)-params.radius:temp(1,2)+params.radius,temp(1,1)-params.radius:temp(1,1)+params.radius)=1;
            all_masks{cur_probe} = mask;
            %flatten
            mask = mask(:);
            mask(nanpxs)=[];        

            dff_probe{cur_file}(:,cur_probe) = nanmean(dff(:,mask==1),2);    
        end    

        %make figures showing the roi location
        if verbose
            figure
            combined_mask = sum(cat(3,all_masks{:}),3);
            ref_img = SpatialBin(opts.cropped_alligned_img,opts.spatial_bin_factor);
            imshowpair(ref_img,combined_mask)   
            [~,fn] = fileparts(dff_fn);
            set(gcf,'position',[681   260   783   460])
            title(fn,'Interpreter','none')
        end
    end %end of imaging data loop
else
   dff_probe = []; 
end

%compile spikeing data
if ~isempty(spike_opts_list)    
    st_probe = cell(1,numel(spike_opts_list));
    n_neurons = [];

    for cur_file = 1:numel(spike_opts_list)    
        fprintf('\n\t working on ephys data for file %d of %d',cur_file, numel(spike_opts_list));
        spike_opts_fn = spike_opts_list{cur_file};
        spike_opts = load(spike_opts_fn);
        spike_opts = spike_opts.opts; 

        %if spock
        if ~ispc
            spike_opts.nidaq_path = ConvertToBucketPath(spike_opts.nidaq_path);
        end
        
        %additional curation
        [ap_clusters,mua_clusters] = CreateMasks(spike_opts); 
        
        %parse spikes
        st_probe{cur_file} = [];
        for cur_probe = 1:numel(spike_opts.kilosort_chan_map_names)
            %load the ephys data
            spikes = load([spike_opts.nidaq_path,sprintf('AP_Probe%d.mat',cur_probe)]);  
            spikes.clust_info.vert_depth = spikes.vert_depth; %add the vertical depth to the structure
            %grab top x mm of spike IDs
            if params.mua == 1 %include labeleds muas
               spikes.clust_info = ApplyMasks(spikes.clust_info,ap_clusters{cur_probe}+mua_clusters{cur_probe});               
            else
               spikes.clust_info = ApplyMasks(spikes.clust_info,ap_clusters{cur_probe}); 
            end
%             %Use angled depth along probe
%             spikes.clust_info.d_surface = SpikeDepth(spikes.clust_info,spike_opts); 
%             IDs = spikes.clust_info.id(spikes.clust_info.d_surface>params.depth(1) & spikes.clust_info.d_surface<params.depth(2));  
            %use vertical depth
            IDs = spikes.clust_info.id(spikes.clust_info.vert_depth>params.depth(1) & spikes.clust_info.vert_depth<params.depth(2));  
            st = spikes.clust_info.spike_times(ismember(spikes.clust_info.spike_cluster,IDs)); 
            n_neurons(cur_file,cur_probe) = numel(IDs);

            %bin to the framerate of the imaging
            fileID = fopen([spike_opts.nidaq_path,'CameraFrameFrontEdgeTimes.txt'],'r');
            formatSpec = '%f';   
            im_times = fscanf(fileID,formatSpec);

            %mua per bin. Binning bins 'to the left' (1st bin = everything
            %between the first and second camera frame', which is the same
            %things that the camera does (i.e. frame one is all the
            %exposure between first edge and second edge)
            st = [histcounts(st,im_times),0]; %since k=nedges-1
               
%             %remove slow fluctuations (This makes no sense, don't do
%             this. only keeping for th efiltering code
%             d = designfilt('highpassiir','FilterOrder',4,'PassbandFrequency',0.02,'PassbandRipple',...
%                 0.1,'Samplerate',floor(nanmedian(diff(im_times))*1000),'DesignMethod','cheby1');
%             st = filtfilt(d,st)+nanmean(st);                
            
            %temporally bin
            if params.bindata ==1  
                st = st(1:2:end)+st(2:2:end);
            end           
            
            %pad with single sample at end to match for ease
            st_probe{cur_file}(:,cur_probe) = st';            
        end
    end %file loop
else
    st_probe = [];
end

%collect figure handles
fig_handles = get(groot, 'Children');

end %function end




















