%outline
%load the cortical traces for each probe
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifTriggeredRaster';
if ~exist(savedir,'dir'); mkdir(savedir); end
figtitle = 'Mouse 331 Recording 1';
ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_11_2021_1dff_combined_processed.mat';
ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_11_2021_probe_coords.mat';
EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
[motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});

%load ephys data
[st_mat,opts,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);

%load the deconvolve dff form the probe insertion sites
[dff,probe_offset] = LoadInsertionSiteDFF(ImgPath,ImgProbeLoc); %probe offset is a few pixel offset in the probe insertion location used for deconvolution to adjust for the fact that the probes are inserted at an angle so the first recorded cell bodies are m/l and a/p shifted from the exact probe insertion site at the surface. 

%load the probe coordinates 
probe_loc = load(ImgProbeLoc); probe_loc = probe_loc.probe_coords;

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%load the motifs 
W_basis = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery\Mouse_basis_motifs.mat','W_basis_Center');
W_basis = W_basis.W_basis_Center; 

%add the motifs panels to the edge of the figure. 

%% make plots
close all;
label = {'RSP/PAR','SS','VIS','M2/PFC'};
win = [-5,15];
for cur_motif = 1:15
    w_max = squeeze(W_basis(:,cur_motif,:));
    [~,idx] = maxk(nanmean(w_max),5);
    w_max = reshape(w_max(:,sort(idx,'ascend')),[68 68 5]);
    
    [trig_dff,trig_st] = ParseByOnset(dff,st_norm,motif_onset,win,cur_motif);
    
    figure('units','normalized','position',[0.2719 0.0657 0.6979 0.8398]); hold on; 
    set(gcf,'color','w')
    t=tiledlayout(5,5,'Padding','normal','TileSpacing','normal');
    title(t,sprintf('Motif %d',cur_motif))
    
    %plot the 5 panels of the motif along the left side
    tileidx = [1,6,11,16,21];
    for i = 1:5
        nexttile(tileidx(i)); hold on; 
        set(gca,'ydir','reverse')
        imagesc(w_max(8:60,4:66,i),[0,prctile(w_max(:),99)]); axis off; axis square;
        colormap(gca,'magma')
        arrayfun(@(n) plot(probe_loc{n}(1,1)+probe_offset{n}(1),probe_loc{n}(1,2)+probe_offset{n}(2),'linewidth',1,'color',col(n,:),'marker','x','markersize',5),1:numel(label),'UniformOutput',0)            
    end
    
    %plot the cortical area activity
    col = getColorPalet(4);
    for i =1:4
       nexttile(i+1); hold on;
       shadedErrorBar(1:size(trig_dff,2),nanmean(trig_dff(i,:,:),3),sem(trig_dff(i,:,:),3),'lineprops',{'color',col(i,:)});   
       ylim([2 5])
       line(get(gca,'xlim'),repmat(nanmean(dff(:,i)),2),'color','k','linestyle','--','linewidth',2); %line denoting avg fr   
       line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset   
       if i == 1
           ylabel('deconvolved FR');
       end
       title(sprintf('%s Probe',label{i}))
    end

    %plot the ephys 
    for i =1:4
        nexttile([4,1]); hold on;
        c = getColorPalet(10); 
        imagesc(nanmean(trig_st{i},3),[0 2.5]); colormap(gca,flipud(gray));
        set(gca,'ydir','reverse'); 
        set(gca,'ylim',[0 size(trig_st{i},1)],'xlim',[0 size(trig_st{i},2)],'YTick','')
        PlotProbeAnatomy(gca, neu_area{i}, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:)); 
        set(gca,'ylim',[0 size(trig_st{i},1)],'xlim',[0 size(trig_st{i},2)],'YTick','')
        line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset  
    end
end
%%











%% motif triggered figures
for  n=1:8  
close all
motif_onset = [onsets{:,n}];
temp = arrayfun(@(x) dff_roi(x-5:x+35,:),motif_onset,'UniformOutput',0);
temp = cat(3,temp{:});
col = getColorPalet(4);
figure; hold on; 
for i =1:4
   shadedErrorBar(1:size(temp,1),nanmean(temp(:,i,:),3),sem(temp(:,i,:),3),'lineprops',{'color',col(i,:)});
end
ylabel({'mean zscore','in craniotomy'})
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
title(sprintf('deconvolved dff, motif %d',n));

% motif activity
[trace,~] = dffROI(squeeze(w_good(:,n,:))',nanpxs,coords);
figure; hold on; 
for i = 1:4
    plot(cat(1,zeros(5,1),trace(:,i)),'color',col(i,:),'linewidth',2);
end
ylabel({'motif activity','in craniotomy'})
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
title(sprintf('Motif dff, motif %d',n));

% plot the average ephys activity in superficial probe
temp = cellfun(@(x) nanmean(x,2), st_mat_super,'UniformOutput',0);
temp = [temp{:}];
temp = temp./std(temp,[],1);
temp = arrayfun(@(x) temp(x-5:x+35,:),motif_onset,'UniformOutput',0);
temp = cat(3,temp{:});
figure; hold on; 
for i =1:4
   shadedErrorBar(1:size(temp,1),nanmean(temp(:,i,:),3),sem(temp(:,i,:),3),'lineprops',{'color',col(i,:)});
end
ylabel({'std-normalized activity'})
title(sprintf('top 600mm, mean spiking, motif %d',n));
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
% plot the motif heatmap organized by depth
temp_mat = cellfun(@(x) zscore(x), st_mat,'UniformOutput',0);
%reorganize by dpeth 
for i = 1:numel(temp_mat)
   [~,idx] = sort(st_depth{i},'ascend');
   temp_mat{i} = temp_mat{i}(:,idx);   
end

figure('position',[108 391 1654 588]); hold on; 
label = {'par','ss','v-thal','m2-pfc'};
for i = 1:numel(temp_mat)
    subplot(1,4,i);
    temp = arrayfun(@(x) temp_mat{i}(x-5:x+35,:),motif_onset,'UniformOutput',0);
    temp = nanmean(cat(3,temp{:}),3);
    imagesc(temp'); colormap gray; colorbar
    title(sprintf('%s \ndepth=%0.2f cm',label{i},max(st_depth{i})/1000))
    sgtitle(sprintf('motif %d',n));
    line([5,5],get(gca,'ylim'),'linewidth',2,'color',col(i,:));
    ylabel('neuron');
end

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\motif_trig_fig';
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('TriggeredActivity_motif%d',n),savedir,0); close all
end



