function Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
%Camden - timeless
%computes loads all data, computes motif onsets and plots the corresponding
%spike rasters. See bottom of function for code to run across all animals
%and recrodings. 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifTriggeredRaster';
if ~exist(savedir,'dir'); mkdir(savedir); end

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
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
    col = getColorPalet(4);
    tileidx = [1,6,11,16,21];
    for i = 1:5
        nexttile(tileidx(i)); hold on; 
        set(gca,'ydir','reverse')
        imagesc(w_max(8:60,4:66,i),[0,prctile(w_max(:),99)]); axis off; axis square;
        colormap(gca,'magma')
        arrayfun(@(n) plot(probe_loc{n}(1,1)+probe_offset{n}(1),probe_loc{n}(1,2)+probe_offset{n}(2),'linewidth',1,'color',col(n,:),'marker','x','markersize',5),1:numel(label),'UniformOutput',0)            
    end
    
    %plot the cortical area activity    
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
        c = getColorPalet(50); 
        imagesc(nanmean(trig_st{i},3),[0 2.5]); colormap(gca,flipud(gray));
        set(gca,'ydir','reverse');         
        PlotProbeAnatomy(gca, neu_area{i}, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:)); 
        set(gca,'ylim',[0 size(trig_st{i},1)],'xlim',[0 size(trig_st{i},2)],'YTick','')
        line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset  
    end
end

[~,temp] = fileparts(ImgPath);
temp = erase(temp,'dff_combined_processed');
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('TrigRaster_rec%s',temp),savedir,0); close all

end %function


% %Code to run on all animals
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_11_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_11_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
% 
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_12_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_12_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
% 
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
% 
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_08_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_08_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
% 
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse334_RestingState_NP_06_09_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_09_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
% 
% ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse334_RestingState_NP_06_10_2021_1dff_combined_processed.mat';
% ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_10_2021_probe_coords.mat';
% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% Plot_MotifTriggeredRasters(ImgPath,ImgProbeLoc,EphysPath,motif_fits)
