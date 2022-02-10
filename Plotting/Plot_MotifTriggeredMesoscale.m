function Plot_MotifTriggeredMesoscale(cur_rec)
%Camden - timeless
%computes loads all data, computes motif onsets and plots the corresponding
%spike rasters. See bottom of function for code to run across all animals
%and recrodings. 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifTriggeredMesoscale';
if ~exist(savedir,'dir'); mkdir(savedir); end

[rectitle,ImgPath,ImgProbeLoc,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%load the imaging 
load(ImgPath,'data_norm','nanpxs');

%% make plots
win = [-1 15];
winval = win(1):win(2);
for cur_motif = 1:15
    %parse by onset 
    trig_dff = ParseByOnset(data_norm',[],motif_onset,win,cur_motif);
    %get average activity
    dff = conditionDffMat(nanmean(trig_dff,3)',nanpxs);
    figure('units','normalized','position',[0 0 1 1]); hold on; 
    [num_rows, num_col]=numSubplot(size(dff,3),0.5);
    cval = [prctile(dff(:),1) prctile(dff(:),99.9)];
    for j = 1:size(dff,3)
       subplot(num_rows,num_col,j); hold on;        
       imagesc(dff(:,:,j),cval); colormap magma;
       axis square; set(gca,'ydir','reverse');
       set(gca,'xlim',[0 68],'ylim',[0 68]); axis off;
       title(sprintf('%d',winval(j)),'fontweight','normal');
       sgtitle(sprintf('motif %d',cur_motif),'fontweight','normal');
       if j==size(dff,3); colorbar; end
    end
    %also plot the variance in residuals 
    dff = conditionDffMat(nanvar(trig_dff-nanmean(trig_dff,3),[],3)',nanpxs);
    figure('units','normalized','position',[0 0 1 1]); hold on; 
    [num_rows, num_col]=numSubplot(size(dff,3),0.5);
    cval = [prctile(dff(:),1) prctile(dff(:),99.9)];    
    for j = 1:size(dff,3)
       subplot(num_rows,num_col,j); hold on; 
       imagesc(dff(:,:,j),cval); colormap magma;
       axis square; set(gca,'ydir','reverse');
       set(gca,'xlim',[0 68],'ylim',[0 68]); axis off;
       title(sprintf('%d',winval(j)),'fontweight','normal');
       sgtitle(sprintf('Variance in Residuals motif %d',cur_motif),'fontweight','normal');
       if j==size(dff,3); colorbar; end
    end  
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('TrigMeso_rec%s_motif%d',rectitle,cur_motif),savedir,0); close all
end


end %function

