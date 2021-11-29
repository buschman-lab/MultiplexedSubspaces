function PlotEstimatedFRDistribution(data,nanpxs,mouse)
%Camden MacDowell - timeless
%Plots the firing rate distribution across different areas
%Data and nanpxs are the preprocessed and deconvolved signal

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Deconvolution\FRDistributionGutchecks';
if ~exist(savedir,'dir'); mkdir(savedir); end
%load fit to 600uM 
fp = fig_params_deconvolutionpaper;

%basic histogram
figure; hold on; 
histogram(data(:));
yvals = get(gca,'ylim');
plot([0 0],yvals,'linestyle',':','color','k','linewidth',2)
xlabel('zcore');
fp.FormatAxes(gca); 
title(sprintf('Distribution all pixels\n Percent negative %.2g\n',100*sum(data(:)<0)/numel(data)));
fp.FigureSizing(gcf,[3 3 6 6],[])


%average histogram per pixel
edges = [-10:0.1:10];
c = NaN(size(data,2),numel(edges)-1);
for i = 1:size(data,2)
   c(i,:) = histcounts(data(:,i),edges);  
end
figure; hold on; 
shadedErrorBar(edges(1:end-1),nanmean(c),std(c),'lineprops',{'color','k','linewidth',2});
yvals = get(gca,'ylim');
plot([0 0],yvals,'linestyle',':','color','k','linewidth',2)
xlim([-2 6])
temp = arrayfun(@(n) sum(data(:,n)<0)/numel(data(:,n)),1:size(data,2));
xlabel('zcore');
fp.FormatAxes(gca); 
title(sprintf('Distribution (mean across pixels)\n mean percent negative %.2g',100*nanmean(temp)));
fp.FigureSizing(gcf,[3 3 6 6],[])

%spatial distribution of negative firing rates
temp = arrayfun(@(n) sum(data(:,n)<0)/numel(data(:,n)),1:size(data,2));
temp = 100*conditionDffMat(temp,nanpxs);
figure; hold on; 
imagesc(temp,[-10,25]); colormap magma; colorbar;
axis off; set(gca,'ydir','reverse')
fp.FormatAxes(gca); 
fp.FigureSizing(gcf,[3 3 6 6],[])
title(sprintf('Spatial distribution of %% \nnegative firing rates'),'fontsize',10);

saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('FR_negative_Gutchecks_rec%s',mouse),savedir,0); close all



end

