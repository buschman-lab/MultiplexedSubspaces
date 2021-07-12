function FitKernel(dff_fn,spike_opts_fn,params)
%Camden MacDowell - timeless

%include MUA
if nargin <3
    params.mua = 0; %1= use both 'good' and 'mua' units. 0 = just 'good'
    params.depth = 800; %depth from surface of probe
    params.radius = 2; %pixel radius around probe tip
    params.offset = { %moves center of DFF circle for probes at angles
end

%load the imaging data
data = load(dff_fn);
dff = data.dff; 
[x,y,z] = size(dff); 
[dff,nanpxs] = conditionDffMat(dff); 
coords = data.probe_coords; 
opts = data.opts; 
clear data

%load ephys opts
spike_opts = load(spike_opts_fn);
spike_opts = spike_opts.opts; 

r = params.radius; %pixel radius around probe tip
spike_d =params.depth; %mm of spike depth

%parse the radius around the probes
dff_probe = NaN(z,numel(coords));
for cur_probe = 1:numel(coords)
    %create mask of radius around probe tip
    temp = coords{cur_probe};
    mask = zeros(x,y);  
    mask(temp(1,2)-r:temp(1,2)+r,temp(1,1)-r:temp(1,1)+r)=1;
    
    %flatten
    mask(nanpxs)=[]; 
    
    dff_probe(:,cur_probe) = nanmean(dff(:,mask==1),2);    
end

%loop through probes
for cur_probe = 1:numel(coords)
   spikes = load([spike_opts.nidaq_path,sprintf('AP_Probe%d.mat',cur_probe)]);   
   
   %grab top x mm of spike IDs
   if params.mua == 1 %include labeleds muas
       idx = 1:numel(spikes.clust_info.depth);       
   else
       idx = sum(ismember(spikes.clust_info.group,'good '),2)==5;       
   end
   IDs = spikes.clust_info.id(spikes.clust_info.depth(idx)<1000);     
   st = spikes.clust_info.spike_times(ismember(spikes.clust_info.spike_cluster,IDs)); 
   
   %bin to the framerate of the imaging
   fileID = fopen([spike_opts.nidaq_path,'CameraFrameFrontEdgeTimes.txt'],'r');
   formatSpec = '%f';   
   im_times = fscanf(fileID,formatSpec);
   
   %mua per bin
   st = histcounts (st,im_times);
   
   %std normalize with a moving window to match dff.
   st = st./movstd(st,opts.method_window*opts.fps);  
   %pad with single sample at end to match for ease
   st = [st,0]';
   
   %look at the cross correlation
   num_iter =1000;  
   binsize = opts.fps*5;
   xcorrkernel = NaN(num_iter,binsize*2+1);   
   dur=(opts.fps*10*60);       
   for i = 1:num_iter   
      idx = randi([dur+1,numel(st)-dur]);
      [xcorrkernel(i,:),t] = xcorr(st(idx:idx+dur),dff_probe(idx:idx+dur,cur_probe)',binsize,'normalized');
   end
   figure; hold on; 
   plot(t/30,nanmean(xcorrkernel,1))
   xval = get(gca,'xlim');
   yval = get(gca,'ylim');
   line([0 0],[-.1 .1],'color','k','linestyle',':')
   title(sprintf('Probe %d Xcorr',cur_probe));
   xlabel('Widefield Lag (s)')
   ylabel('Correlation');
   xlim([xval]); ylim([yval]);     
   
   %fit a regression kernel where you predict fr at each timepoint using
   %imaging data from subsequent 2 second   
   
   binsize = opts.fps*2;
   kernel=[];
   for j = 1:1000
       str = randi([2*binsize,numel(st)-(2*opts.fps*60)]);
       count=1;
       predictors = [];
       response = [];   
       for i = str:str+opts.fps*60
           predictors(count,:)= dff_probe(i-binsize:i+binsize,cur_probe)';
           response(count) = st(i);    
           count= count+1;
       end
       kernel(j,:) = glmfit(predictors,response);
   end
   %remove intercept
   kernel = kernel(:,2:end);
   %maximum normalized
   kernel_norm = kernel./max(abs(kernel),[],2);
   kernel_mean = nanmean(kernel_norm,1);       
   kernel_mean = kernel_mean./sum(kernel_mean.^2);
   figure; hold on;
   t=[-1*binsize:binsize]/30;   
   plot(t,kernel_mean,'color','k','linewidth',2);   
   yval = get(gca,'ylim');   
   line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
   xlim([-0.75 0.75]); ylim([yval]);
   title(sprintf('Probe %d Kernel',cur_probe));
   xlabel('flourescence Lag (s)');     
end %probe loop

%save off figures
handles = get(groot, 'Children'); 
file_list={};
for i = 1:numel(handles)
    file_list{i} = [pwd filesep ,sprintf('%d_kernel.pdf',i)];
    print(handles(i),file_list{i},'-dpdf','-bestfit');
end
append_pdfs(sprintf('pca_kernels_muas%d_d%d_r%d.pdf',params.mua,params.depth,params.radius),file_list{:});

% %delete original figures
cellfun(@(x) delete(x),file_list)


end %function



%take top X mm spikes (and maybe MUA)

%load imaging frame times 

%bin and std normalize spikes

% 
%    binsize = opts.fps*2;
%    num_iter =1000;      
%    predictors = NaN(num_iter,binsize*2);
%    response = NaN(num_iter,1);
%    for i = 1:num_iter
%        idx = randi([binsize+1,numel(st)-binsize]);
%        predictors(i,:) = dff_probe(idx-binsize:idx+binsize-1,cur_probe)';
% %        response(i) = st(i);
%        response(i) = dff_probe(i,cur_probe)';
% %        predictors(i,:) = st(idx-binsize:idx+binsize-1)';
%    end
%    mdl= fitlm(predictors,response);
%    kernel = mdl.Coefficients{2:end,1};
%    figure; plot(kernel);

%     kernel_t = [-0.5,0.5];
%     kernel_frames = floor(kernel_t(1)*sample_rate):ceil(kernel_t(2)*sample_rate);  
%     zs = [false,false];
%     cvfold = 5;
%     return_constant = true;
%     use_constant = true;
%     lambda = 0;


%fit kernel


%compare across locations
