function FitKernel(dff_fn,spike_opts_fn,params)
%Camden MacDowell - timeless

%include MUA
if nargin <3
    params.mua = 0; %1= use both 'good' and 'mua' units. 0 = just 'good'
    params.depth = 800; %depth from surface of probe
    params.radius = 2; %pixel radius around probe tip
    params.offset = {[1,0],[0,-1],[-1,1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
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

%loop through probes
for cur_probe = 2%1:numel(coords)
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
   st = histcounts (st,im_times);
   
   %std normalize with a moving window to match dff.
   st = st./movstd(st,opts.method_window*opts.fps);  
   %pad with single sample at end to match for ease
   st = [st,0]';
   
   %img --> ap xcorr
   maxlag= opts.fps*2.5;
   [xcorrkernel,t] = xcorr(st-nanmean(st),dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
   st_auto = xcorr(st-nanmean(st),maxlag,'normalized');
   img_auto = xcorr(dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
   figure; hold on; 
   plot(t/30,xcorrkernel,'k','linewidth',2); 
   xval = get(gca,'xlim');
   yval = get(gca,'ylim');
   line([0 0],[-.1 .1],'color','k','linestyle',':')
   title(sprintf('Probe %d Xcorr',cur_probe));
   xlabel('Widefield Lag (s)')
   ylabel('Correlation');
   xlim([xval]); ylim([yval]);    
   
   %img-->img and ap-->ap autocorr
   figure; hold on; 
   p1=plot(t/30,st_auto,'r','linewidth',2,'linestyle',':');
   p2=plot(t/30,img_auto,'b','linewidth',2,'linestyle','--');
   xval = get(gca,'xlim');
   yval = get(gca,'ylim');
   line([0 0],[-.1 .1],'color','k','linestyle',':')
   title(sprintf('Probe %d Autocorr',cur_probe));
   xlabel('Lag (s)')
   legend([p1,p2],'spikes','img');
   ylabel('Correlation');
   xlim([xval]); ylim([yval]);   
   
   %img --> ap kernel | predict fr at each timepoint using 1s imaging data    
   %individual points
   win = opts.fps*2;
   if mod(win,2)==0; win = win+1;  end %odd so can appropriately center response within window
   predictors = createRollingWindow(dff_probe(:,cur_probe), win);
   response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
%    predictors = createRollingWindow(st, win);
%    response =  dff_probe(ceil(win/2):end-floor(win/2),cur_probe); % get the middle timepoitn in window     
   kernel = NaN(numel(win),1);
   for i = 1:size(predictors,2)
      temp = glmfit(predictors(:,i),response); %mean centering not needed here since you don't care about intercept
      kernel(i,1) = temp(2);
   end
%    kernel = kernel(2:end); %remove intercept. 
      
   figure; hold on;
   t=[-1*floor(win/2):floor(win/2)]/30;   
   plot(t,kernel,'color','k','linewidth',2);   
   yval = get(gca,'ylim');   
   line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
   xlim([-0.75 0.75]); ylim([yval]);
   title(sprintf('Probe %d Kernel Single Points',cur_probe));
   xlabel('flourescence Lag (s)');    
   
      %kernel full window
   if mod(win,2)==0; win = win+1;  end %odd so can appropriately center response within window
   predictors = createRollingWindow(dff_probe(:,cur_probe), win);
   response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
   kernel= glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
   kernel = kernel(2:end); %remove intercept. 
      
   figure; hold on;
   t=[-1*floor(win/2):floor(win/2)]/30;   
   plot(t,kernel,'color','k','linewidth',2);   
   yval = get(gca,'ylim');   
   line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
   xlim([-0.75 0.75]); ylim([yval]);
   title(sprintf('Probe %d Kernel Full Window',cur_probe));
   xlabel('flourescence Lag (s)'); 
    
   %test how well it works by looks at the autocorrelation before and after
   dff_conv = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid');
   figure; hold on; 
   subplot(2,3,1)
   histogram(dff_probe(:,cur_probe));
   subplot(2,3,2)
   histogram(dff_conv(:));
   subplot(2,3,3)
   histogram(st(:));
   subplot(2,3,[4,6]); hold on;
   plot(dff_probe(1:opts.fps*30,cur_probe),'color',[0.5 0.5 0.5 0.5],'linewidth',2)
   plot(dff_conv(1:opts.fps*30),'color','k','linewidth',2)
   yyaxis right
   plot(st(1:opts.fps*30),'color','r','linewidth',2,'linestyle',':')
   
   %deconvolved img --> ap xcorr and deconv img --> deconv --img
   [xcorrkernel_deconv,t] = xcorr(st-nanmean(st),dff_conv'-nanmean(dff_conv),maxlag,'normalized');   
   xcorrkernel = xcorr(st-nanmean(st),dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
   deconv_auto = xcorr(dff_conv'-nanmean(dff_conv),maxlag,'normalized');
   img_auto = xcorr(dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
   figure; hold on; 
   p0=plot(t/30,xcorrkernel_deconv,'k','linewidth',2);    
   p1=plot(t/30,xcorrkernel,'k','linewidth',2,'linestyle',':');    
   xval = get(gca,'xlim');
   yval = get(gca,'ylim');
   line([0 0],[-.1 .1],'color','k','linestyle',':')   
   xlabel(' Lag (s)')
   ylabel('Correlation');
   xlim([xval]); ylim([yval]);   
   yyaxis right
   p2=plot(t/30,deconv_auto,'g','linewidth',2,'linestyle','--');
   p3=plot(t/30,st_auto,'r','linewidth',2,'linestyle',':');
   p4=plot(t/30,img_auto,'b','linewidth',2,'linestyle','--');
   legend([p0,p1,p2,p3,p4],'deconv-->st','img-->st','deconv-->deconv','st-->st','img-->img');
   set(gca,'ycolor','k'); ylabel('Autocorrelation');
   
   %include in title the quality of deconv fit to fr
   title(sprintf('Probe %d Xcorr Deconv RMSE=%0.2g Img RMSE=%0.2g',cur_probe,sqrt(mse(st-dff_conv)),sqrt(mse(st-dff_probe(:,cur_probe)))));
      
   %ap --> ap autokernel img --> img autokernel 
   predictors = createRollingWindow(st, win);
%    predictors(:,ceil(win/2)-5:ceil(win/2)+5)=0;
   response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window   
   kernel = glmfit(predictors,response);
   kernel = kernel(2:end); %remove intercept. 
      
   figure; hold on;
   t=[-1*floor(win/2):floor(win/2)]/30;   
   plot(t,kernel,'color','k','linewidth',2);   
   yval = get(gca,'ylim');   
   line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
   xlim([-0.75 0.75]); ylim([yval]);
   title(sprintf('ST GLM d%d r%d MUA%d',params.depth,params.radius,params.mua));
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

% add in self prediction for both imaging and neural
% remake so using the entire data
% add probe adjustments


%%old
% num_iter =1000;  
% binsize = opts.fps*2.5;
% xcorrkernel = NaN(num_iter,binsize*2+1);   
% dur=(opts.fps*10*60);       
% for i = 1:num_iter   
%   idx = randi([dur+1,numel(st)-dur]);
%   [xcorrkernel(i,:),t] = xcorr(st(idx:idx+dur),dff_probe(idx:idx+dur,cur_probe)',binsize,'normalized');
% end


%    binsize = opts.fps*2;
%    kernel=[];
%    for j = 1:1000
%        str = randi([2*binsize,numel(st)-(2*opts.fps*60)]);
%        count=1;
%        predictors = [];
%        response = [];   
%        for i = str:str+opts.fps*60
%            predictors(count,:)= dff_probe(i-binsize:i+binsize,cur_probe)';
%            response(count) = st(i);    
%            count= count+1;
%        end
%        kernel(j,:) = glmfit(predictors,response);
%    end
%    %remove intercept
%    kernel = kernel(:,2:end);
%    %maximum normalized
%    kernel_norm = kernel./max(abs(kernel),[],2);
%    kernel_mean = nanmean(kernel_norm,1);       
%    kernel_mean = kernel_mean./sum(kernel_mean.^2);