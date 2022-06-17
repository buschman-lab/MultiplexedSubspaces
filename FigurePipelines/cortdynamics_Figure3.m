%Figure 3
fp = fig_params_cortdynamics;
%load data
data = LoadCorticalNetworks(0);
dataout = LoadCorticalNetworks(1);

%% Plot example Cortical Networks

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ExampleCorticalNetworks';
if ~exist(savedir,'dir'); mkdir(savedir); end
binarized = 1 ; %def = 0
threshold = 0.25; %def = 0.25
for cur_rec = 1 %make sure results are consistent across recordings
    cur_motif = 8; cur_area = 2; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,98) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all
    
    cur_motif = 11; cur_area = 7; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,98) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all
    
    cur_motif = 6; cur_area = 4; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,98) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all 
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all
    
    cur_motif = 14; cur_area = 3; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,98) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all

    cur_motif = 5; cur_area = 8; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,97) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all

    cur_motif = 9; cur_area = 8; 
%     Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,1,98) %(data,cur_rec,motif,area,sigflag)
%     saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all
    CorticalHexPlots(data,cur_rec,cur_motif,cur_area,'Threshold',threshold,'Binarize',binarized); 
    title(sprintf('Binarized %d Thresh %0.2f',binarized,threshold),'fontweight','normal');
    saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_barcode_thres%.0f',cur_motif,cur_rec,cur_area,threshold*10),savedir,0); close all
    
end


%% Plot the overlap (histogram) of subspace networks 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\CorticalNetworkSummaryStats';
Plot_OverlapAcrossDimensions(data,1); %do spatial correlation
Plot_OverlapAcrossDimensions(data,0); %do spatial overlap
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('WithinSubspaceOverlap'),savedir,0); close all
    
%% Plot how many network individual neural pixels participated in
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\CorticalNetworkSummaryStats';
Plot_NumberOfDimension(data); %do spatial correlation
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('DimensionsPerPixel'),savedir,0); close all


%% Plot the overall strength and the size of the significance-thresholded network across dimensions
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\CorticalNetworkSummaryStats';
if ~exist(savedir); mkdir(savedir); end
Plot_CortNetworksSummaryStats(data)
saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},'SumStats',savedir,0); close all


%% Plot the overlap between dimensions across motifs
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\CorticalNetworkSummaryStats';
if ~exist(savedir); mkdir(savedir); end
Plot_MotifCorticalNetworkOverlap(data)
saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},'XMotifOverlap',savedir,0); close all





%% example across animals
%vis and motif 5
cur_motif = 5; cur_rec = 3; cur_area = 8; 
Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,0,98) %(data,cur_rec,motif,area,sigflag)
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d',cur_motif,cur_rec,cur_area),savedir,0); close all
cur_motif = 5; cur_rec = 5; cur_area = 8; 
Plot_ExampleCorticalNetworks(data,cur_rec,cur_motif,cur_area,0,98) %(data,cur_rec,motif,area,sigflag)
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('motif%d_rec%d_area%d_sig',cur_motif,cur_rec,cur_area),savedir,0); close all


%Show that they are the same across animals

%% Plot the clusters
% num_clust = Plot_CorticalNetworkClusters(data,'wRegionxMotif');
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\CorticalNetworkSummaryStats';
num_clust = Plot_CorticalNetworkClusters(data,'overall',1);
saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},sprintf('SumStatsSimilarityMatrix',cur_motif,cur_rec,cur_area),savedir,0); close all
set(gca,'clim',[.2,1])
%% Plot the size of the network with increasing dimensions
% Plot_CorticalNetworkSize(data);



%% Analyze data
[~,arearho,arearho_ordered,areaorder,arearho_avg,inoutrho,inoutrho_avg,xarearho,xarearho_ordered,xarearho_avg] = AnalyzeCorticalNetworks(data,dataout);

%basically for each spatial distribution I want to plot the max similarity
%with other motifs. But then show that they are reordered. 

%% Plot the ordered and non_ordered correlation for an area
close all
area_label = data{1}(1).area_label;
x = fisherZ(arearho-arearho_avg);
x = squeeze(nanmean(x,1));
xx = fisherZ(arearho_ordered-arearho_avg);
xx = squeeze(nanmean(xx,1));

col = fp.c_area;
ndim=10; 
for i = 1:size(col,1)
    figure; hold on;
    y = squeeze(x(i,:,:));
%     shadedErrorBar(1:ndim,nanmean(y),sem(y,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});    
    yy = squeeze(xx(i,:,:));
%     shadedErrorBar(1:ndim,nanmean(yy),sem(yy,1),'lineprops',{'color',[col(i,:),0.75],'linestyle',':','linewidth',2});
    arrayfun(@(n) plot(1:size(y,2), y(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.75,'markersize',4), 1:size(y,1));
    plot(1:size(y,2), nanmean(y),'linestyle','-','marker','none','color',col(i,:),'linewidth',2.5,'markersize',5,'MarkerFaceColor','none');

    arrayfun(@(n) plot(1:size(yy,2), yy(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.75,'markersize',4), 1:size(yy,1));
    plot(1:size(yy,2), nanmean(yy),'linestyle',':','marker','none','color',col(i,:),'linewidth',2.5,'markersize',5,'MarkerFaceColor','none');
    xlim([1 ndim]);
    xlabel('subspace dimension');
    ylabel({'Spatial Rho','(global average subtracted)'});
    title(sprintf('Similarity Between \nMotifs |%s ',area_label{i}),'Fontweight','normal')
    legend('best match','ordered','location','bestoutside')
    fp.FormatAxes(gca);  box on; grid on
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 15 10])
end

%% Same but with In/Out as reference
close all
area_label = data{1}(1).area_label;
x = fisherZ(arearho-arearho_avg);
x = squeeze(nanmean(x,1));
xx = fisherZ(arearho_ordered-arearho_avg);
xx = squeeze(nanmean(xx,1));
xxx = fisherZ(inoutrho-inoutrho_avg);
xxx = squeeze(nanmean(xxx,1));
% ii = squeeze(nanmean(areaorder,1));
col = fp.c_area;
ndim=10; 
for i = 1:size(col,1)
    figure; hold on;
    y = squeeze(x(i,:,:));
    shadedErrorBar(1:ndim,nanmean(y),sem(y,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});    
    yy = squeeze(xx(i,:,:));
    shadedErrorBar(1:ndim,nanmean(yy),sem(yy,1),'lineprops',{'color',[col(i,:),0.75],'linestyle',':','linewidth',2});
    yyy = squeeze(xxx(i,:,:));
    shadedErrorBar(1:ndim,nanmean(yyy),sem(yyy,1),'lineprops',{'color',[0.25, 0.25, 0.25, 0.75],'linestyle',':','linewidth',2});            
    xlim([1 ndim]);
    xlabel('subspace dimension');
    ylabel({'Spatial Rho','(global average subtracted)'});
    title(sprintf('Similarity Between \nMotifs |%s ',area_label{i}),'Fontweight','normal')
    legend('best match','ordered','location','bestoutside')
%     yyaxis right
%     tt = squeeze(ii(i,:,:));
%     shadedErrorBar(1:ndim,nanmean(tt),sem(tt,1),'lineprops',{'color',[0.25 0.25 0.25 0.5],'linestyle','-','linewidth',2});        
    fp.FormatAxes(gca);  box on; grid on
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 15 10])
end

%% Plot the correlation across areas within a motif
close all
area_label = data{1}(1).area_label;
x = fisherZ(xarearho-xarearho_avg);
x = squeeze(nanmean(x,1));
xx = fisherZ(xarearho_ordered-xarearho_avg);
xx = squeeze(nanmean(xx,1));
col = fp.c_area;
ndim=10; 
for i = 1:size(col,1)
    figure; hold on;
    y = squeeze(x(i,:,:));
    shadedErrorBar(1:ndim,nanmean(y),sem(y,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});    
    yy = squeeze(xx(i,:,:));
    shadedErrorBar(1:ndim,nanmean(yy),sem(yy,1),'lineprops',{'color',[col(i,:),0.75],'linestyle',':','linewidth',2});
    xlim([1 ndim]);
    xlabel('subspace dimension');
    ylabel({'Spatial Rho','(global average subtracted)'});
    title(sprintf('Similarity Between \nAreas |%s ',area_label{i}),'Fontweight','normal')
    legend('best match','ordered','location','bestoutside')
    fp.FormatAxes(gca);  box on; grid on
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 15 10])
end


%%

% % Other random plots
% % %% Plot the average autocorrelations
% % %average across recordings
% % area_label = data{1}(1).area_label;
% % x = fisherZ(autorho);
% % x = squeeze(nanmean(x,1));
% % col = fp.c_area;
% % ndim=10;
% % figure; hold on; 
% % for i = 1:size(col,1)
% %     temp = squeeze(x(i,:,:));
% %     shadedErrorBar(1:ndim-1,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(i,:),0.25],'linewidth',2});
% % end
% % legend(area_label,'location','bestoutside')
% % xlim([1 ndim-1]);
% % xlabel('subspace dimension');
% % ylabel('Spatial Rho_z');
% % title(sprintf('Autocorrelation of cortical networks'),'Fontweight','normal')
% % fp.FormatAxes(gca);  box on; grid on
% % fp.FigureSizing(gcf,[3 2 4 4],[10 10 15 10])