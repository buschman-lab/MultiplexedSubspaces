function Plot_BehavioralStateDecoder(data)

if nargin <1 
   data = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\LeaveOneOutClassification\LeaveOneOut.mat');   
end

fp = fig_params;

%get the contribution of each motif (full - leave-one-out)
contrib = (data.auc_full - data.auc_leave)*100;

%%
close all
figure('Position',[95   137   972   398]); hold on; 
yvals = [-2, 0, 8];
yvals_norm = 1-(yvals-min(yvals))/(max(yvals)-min(yvals));
imagesc(contrib,[yvals(1),yvals(end)])
cmap = customcolormap(yvals_norm,{'#00d5ff','#000000','#ff00bb'});
colormap(cmap)
c = colorbar;
set(c,'YTick',yvals);
ylabel(c,{'Normalized Motif';'Intensity'},'FontSize',16,'Fontweight','normal','FontName','Arial');
set(c,'units','centimeters','position',[18.25 6 0.5 2])
set(gca,'Units','centimeters','Position',[8 6 10 0.75]); 
set(gca,'YTickLabel',[],'XTickLabel',[]);
xlim([0.5 14.5])
set(gca,'YColor','w','XColor','w')
setFigureDefaults

for x_grid = 0.5:1:14+0.5
    line([x_grid,x_grid],[0.5,1.5],'linewidth',1.5,'color','w')    
end

for i = 1:numel(contrib)
   text([i,i],[2,2],sprintf('%0.2g%%',contrib(i)),'Color',[0.1 0.1 0.1],'FontSize',18,'FontWeight','normal','HorizontalAlignment','left','Rotation',90);
end
%%
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ManuscriptRevisionFigures_currentbio';
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','behavioraldecoder',savedir,1);
close all;
