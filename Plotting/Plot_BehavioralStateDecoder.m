function Plot_BehavioralStateDecoder()
%Camden MacDowell - timeless
%contrib is a motif x replication matrix of contributions of each motif to
%the decoding accuracy on withheld data

fn = GrabFiles('LeaveOneOut_holdout_20_replication_',0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\LeaveOneOutClassification'}); 
data_full = cellfun(@(x) load(x,'auc_full'),fn,'UniformOutput',1);
data_leave = cellfun(@(x) load(x,'auc_leave'),fn,'UniformOutput',1);
data_full = cat(1,data_full(:).auc_full);
data_leave = cat(1,data_leave(:).auc_leave);
contrib_temp = arrayfun(@(n) (data_full(n)-data_leave(n,:))*100,1:numel(data_full),'UniformOutput',0);   
contrib_temp = cat(1,contrib_temp{:});
contrib = nanmean(contrib_temp);
contrib_ci = bootci(100,@nanmean,contrib_temp);

motif_id = 1:14;
[~,idx] = sort(contrib,'ascend');
motif_id = motif_id(idx);
contrib = contrib(idx);
contrib_ci =contrib_ci(:,idx); 

%% Plot a line
figure; hold on;
line([0,0],[0 14],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2)
plot(contrib,1:numel(contrib),'o','color','k','linewidth',2)
errorbar(contrib,1:numel(contrib),[],[],contrib-contrib_ci(1,:),contrib_ci(2,:)-contrib,'color','k','linewidth',1.5,'linestyle','none')
xlabel({'Contibution';'(AUC)'})
ylabel('Basis Motif (Sorted)');
set(gca,'ytick',1:numel(contrib),'yticklabel',motif_id);
setFigureDefaults
set(gca,'Units','centimeters','Position',[2 3 3 9.5]); 
set(gcf,'position',[680   400  722   600]);

%%

%% Plot as a tree
fp = fig_params;
figure('Position',[95   137   972   398]); hold on; 
yvals = [-1, 0, 3];
yvals_norm = 1-(yvals-min(yvals))/(max(yvals)-min(yvals));f
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
% yvals = [-2, 0, 8];
% yvals_norm = 1-(yvals-min(yvals))/(max(yvals)-min(yvals));
% imagesc(contrib,[yvals(1),yvals(end)])
% cmap = customcolormap(yvals_norm,{'#00d5ff','#000000','#ff00bb'});

% yvals = [0, 4.5];
% yvals_norm = 1-(yvals-min(yvals))/(max(yvals)-min(yvals));
% imagesc(contrib,[yvals(1),yvals(end)])
% cmap = customcolormap(yvals_norm,{'#00d5ff','#000000'});%'#00d5ff'