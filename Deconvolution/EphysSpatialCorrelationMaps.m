function EphysSpatialCorrelationMaps()
%Camden MacDowell - timeless
%for each electrode site identifies the best roi to use to deal with
%shuffled angles in probe
load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Deconvolution\restingstate_processed_fn.mat','dff_list','spike_opts_list') 
save_dir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\SupplementalFigures\Deconvolution\';

bindata=0; 
%load ephys data
params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
params.depth = [0 500]; %depth from surface of probe
params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
[~,st] = CompileData_deconvolution([],spike_opts_list,params,0);

for cur_file = 1:numel(dff_list) %loop through files
    spikes = st{cur_file}; 
    
    %load full imaging data
    data = load(dff_list{cur_file},'dff','nanpxs','probe_coords');
    dff = data.dff;
    if bindata ==1
       spikes = spikes(1:2:end,:)+spikes(2:2:end,:);  
       dff = dff(1:2:end,:)+dff(2:2:end,:);             
       dff = conditionDffMat(dff,data.nanpxs,[]);
       dff = DenoisePCA(dff);
       [dff,nanpxs] = conditionDffMat(dff);        
    else 
        nanpxs = data.nanpxs;
    end
    coords = data.probe_coords;
    
    spike_opts = load(spike_opts_list{cur_file},'opts');
    run_name = spike_opts.opts.run_name;

    %get dimesions of the full data
    x = sqrt(size(dff,2)+numel(nanpxs));   
    
    %generate correlation map within a +/-100ms window
    for cur_probe = 1:size(spikes,2)
       rho = arrayfun(@(x) xcorr(spikes(:,cur_probe)-nanmean(spikes(:,cur_probe)),dff(:,x)-nanmean(dff(:,x)),3,'normalized'),1:size(dff,2),'UniformOutput',0);
       [rho_best, lag_best] = cellfun(@max , rho, 'UniformOutput',1);
       lag_best = lag_best-4; %zero lag
       rho_zero = corr(spikes(:,cur_probe)-nanmean(spikes(:,cur_probe)),dff-nanmean(dff));       
       maps = conditionDffMat(cat(1,rho_zero,rho_best,lag_best),nanpxs);
       
       %identify the best radius dimensions around the probe using the zero lag
       rho_roi = NaN(1,3);
       for r = 1:3
          temp = coords{cur_probe};
          temp(1,1)=temp(1,1)+params.offset{cur_probe}(1);
          temp(1,2)=temp(1,2)+params.offset{cur_probe}(2);
          mask = zeros(x,x);  
          mask(temp(1,2)-r:temp(1,2)+r,temp(1,1)-r:temp(1,1)+r)=1;
          mask = mask(:);
          mask(nanpxs)=[];
          rho_roi(r) = nanmean(rho_zero(mask==1));
       end
       
       [~,best_rad] = max(rho_roi);       
       figure; 
       subplot(1,3,1); imagesc(maps(:,:,1)); colorbar; colormap magma; set(gca,'ydir','reverse');   hold on;         
       plot(coords{cur_probe}(1,1),coords{cur_probe}(1,2),'marker','.','markersize',20,'color','c'); %add the probe location
       axis off; title(sprintf('Zero Rho R%d P%d Best Radius=%0.0f',cur_file,cur_probe,best_rad));

       subplot(1,3,2); imagesc(maps(:,:,2)); colorbar; colormap magma; set(gca,'ydir','reverse');   hold on;         
       plot(coords{cur_probe}(1,1),coords{cur_probe}(1,2),'marker','.','markersize',20,'color','c'); %add the probe location
       axis off; title(sprintf('Max xRho R%d P%d',cur_file,cur_probe));           

       subplot(1,3,3); imagesc(maps(:,:,3)); colorbar; colormap magma; set(gca,'ydir','reverse');   hold on;         
       plot(coords{cur_probe}(1,1),coords{cur_probe}(1,2),'marker','.','markersize',20,'color','c'); %add the probe location
       axis off; title(sprintf('xRho Lag R%d P%d',cur_file,cur_probe));                      
       set(gcf,'position',[ 65         559        1747         420])
       
       %save off
       saveas(gcf,[save_dir,sprintf('%d_%s_probe%d_correlation_map.png',cur_file,run_name,cur_probe)]);
       close
       
    end   
end        
     
    
    
end



    