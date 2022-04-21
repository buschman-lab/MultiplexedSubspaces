function statresults = PlotSubspaceDimensionality(data,dataout,area_name,cur_rec)
fp = fig_params_cortdynamics;
thresh = 0.8;

%all
if isempty(area_name) %ccombined across areas
    [rrr_d,~] = LoadVariable(data,'rrr_dim',[],thresh);
    col = fp.c_area; 
    localdim = LoadVariable(data,'svca',[],thresh);
    %average across recordings
    rrr_d = squeeze(nanmean(rrr_d,1));
    localdim = squeeze(nanmean(localdim,1));
    Plot_CorrelateValuesBetweenRecordings(rrr_d',localdim','combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
        'ylabel','Local dimensionality','color_flag',1,'corrtype','kendall','col',col,'addjitter',3)
    plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
    xlim([0 15])
    ylim([0 50])
    title({'IN Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 
    fp.FigureSizing(gcf,[3 2 5 5],[10 10 20 10])

else

[rrr_d,area_label] = LoadVariable(data,'rrr_dim',area_name,thresh);
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
localdim = LoadVariable(data,'svca',area_name,thresh);
Plot_CorrelateValuesBetweenRecordings(rrr_d,localdim,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
    'ylabel','Local dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
xlim([0 max(cat(1,rrr_d(:),localdim(:)))])
ylim([0 max(cat(1,rrr_d(:),localdim(:)))])
title({'IN Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%Also do for a single example recording
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
rrr_d = rrr_d(cur_rec,:);
localdim = localdim(cur_rec,:);
Plot_CorrelateValuesBetweenRecordings(rrr_d,localdim,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
    'ylabel','Local dimensionality','color_flag',0,'corrtype','kendall','col',col,'addjitter',3)
plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
xlim([0 max(cat(1,rrr_d(:),localdim(:)))])
ylim([0 max(cat(1,rrr_d(:),localdim(:)))])
title({'IN Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%do the same thing for data out
[rrr_d,area_label] = LoadVariable(dataout,'rrr_dim',area_name,thresh);
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
localdim = LoadVariable(dataout,'svca',area_name,thresh);
Plot_CorrelateValuesBetweenRecordings(rrr_d,localdim,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
    'ylabel','Local dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
xlim([0 max(cat(1,rrr_d(:),localdim(:)))])
ylim([0 max(cat(1,rrr_d(:),localdim(:)))])
title({'OUT Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%do the same thing for data out
[rrr_d,area_label] = LoadVariable(dataout,'rrr_dim',area_name,thresh);
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
localdim = LoadVariable(dataout,'svca',area_name,thresh);
Plot_CorrelateValuesBetweenRecordings(rrr_d,localdim,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
    'ylabel','Local dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
xlim([0 max(cat(1,rrr_d(:),localdim(:)))])
ylim([0 max(cat(1,rrr_d(:),localdim(:)))])
title({'OUT Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%% plot the ratio of local to global across brain areas (paired)
[rrr_d,area_label] = LoadVariable(data,'rrr_dim',[],thresh);
localdim = LoadVariable(data,'svca',[],thresh);
col = fp.c_area; 
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);

%flatten per area
localdim = reshape(localdim,size(localdim,1)*size(localdim,2),size(localdim,3));
rrr_d = reshape(rrr_d,size(rrr_d,1)*size(rrr_d,2),size(rrr_d,3));
temp = localdim./rrr_d;
figure; hold on;
vp = CompareViolins(temp',fp,'col',col,'connectline',[0.25 0.25 0.25 0.50],'plotspread',0,'divfactor',3);
fp.FigureSizing(gcf,[3 2 6 4],[10 10 20 10]); 
set(gca,'XTickLabel',area_label,'XTickLabelRotation',45)
fp.FormatAxes(gca); box on
ylabel('Local D / Sub D');

%Do for each area across recordings
[rrr_d,area_label] = LoadVariable(data,'rrr_dim',[],thresh);
localdim = LoadVariable(data,'svca',[],thresh);
for i = 1:8
    col = fp.c_area; col = col(strcmp(area_label,area_label{i}),:);
    rrr_dtemp = rrr_d(:,:,strcmp(area_label,area_label{i}));
    localdimtemp = localdim(:,:,strcmp(area_label,area_label{i}));
    Plot_CorrelateValuesBetweenRecordings(rrr_dtemp,localdimtemp,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
        'ylabel','Local dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
    plot([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],[0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],'linestyle',':','color','k')
    xlim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    ylim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    title(sprintf('Out | %s',area_label{i}),'fontweight','normal') 
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])
end

%In versus out
[rrr_d,area_label] = LoadVariable(dataout,'rrr_dim',[],thresh);
localdim = LoadVariable(dataout,'svca',[],thresh);
for i = 1:8
    col = fp.c_area; col = col(strcmp(area_label,area_label{i}),:);
    rrr_dtemp = rrr_d(:,:,strcmp(area_label,area_label{i}));
    localdimtemp = localdim(:,:,strcmp(area_label,area_label{i}));
    Plot_CorrelateValuesBetweenRecordings(rrr_dtemp,localdimtemp,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},...
        'ylabel','Local dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
    plot([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],[0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],'linestyle',':','color','k')
    xlim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    ylim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    title(sprintf('Out | %s',area_label{i}),'fontweight','normal') 
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])
end

%% correlation between in and out subspace dimensionality per area
%Do for each area across recordings
[rrr_d,area_label] = LoadVariable(data,'rrr_dim',[],thresh);
[localdim,area_label] = LoadVariable(dataout,'rrr_dim',[],thresh);
for i = 1:8
    col = fp.c_area; col = col(strcmp(area_label,area_label{i}),:);
    rrr_dtemp = rrr_d(:,:,strcmp(area_label,area_label{i}));
    localdimtemp = localdim(:,:,strcmp(area_label,area_label{i}));
    Plot_CorrelateValuesBetweenRecordings(rrr_dtemp,localdimtemp,'combo',fp,'right','xlabel',{'In Dimensionality'},...
        'ylabel','Out dimensionality','color_flag',2,'corrtype','kendall','col',col,'addjitter',3)
    plot([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],[0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))],'linestyle',':','color','k')
    xlim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    ylim([0 max(cat(1,rrr_dtemp(:),localdimtemp(:)))])
    title(sprintf('In versus out | %s',area_label{i}),'fontweight','normal') 
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])
end

%% averaged across recordings
x = squeeze(nanmean(rrr_d,1));
y = squeeze(nanmean(localdim,1));
xy = abs(x-y);
[~,idxsort] = sort(nanmedian(xy,1),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(xy(:,idxsort),'Notch','off')
a = get(get(gca,'children'),'children');
t = get(a,'tag'); 
idx = find(strcmp(t,'Box')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', col(i,:),'LineWidth',1.5); 
end
idx = find(strcmp(t,'Median')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', 'k','LineWidth',1.5); 
end
xlabel('Area');
set(gca,'xticklabel',area_label(idxsort),'XTickLabelRotation',45)
ylabel({'\Delta Dimensionality','(In - Out)'});
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 10 10])

[p,~,stats] = kruskalwallis(xy(:,idxsort),[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
statresults = multcompare(stats,'Display','off');

%Compare dimensionality across regions
[D,area_label] = LoadVariable(data,'rrr_dim',[]);
%do it relative to the local dimensionality
[L,area_label] = LoadVariable(data,'svca',[],0.8);

d = D./L;
%average across recordings
d = squeeze(nanmean(d,1))';

%reorganize by decreasing median
[~,idxsort] = sort(nanmedian(d,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(d(idxsort,:)','Notch','off')
a = get(get(gca,'children'),'children');
t = get(a,'tag'); 
idx = find(strcmp(t,'Box')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', col(i,:),'LineWidth',1.5); 
end
idx = find(strcmp(t,'Median')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', 'k','LineWidth',1.5); 
end
xlabel('Area');
set(gca,'xticklabel',area_label(idxsort),'XTickLabelRotation',45)
ylabel({'Relative dimensionality','(subspace/local)'});
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])

[p,~,stats] = kruskalwallis(d(idxsort,:)',[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
figure; statresults = multcompare(stats,'Display','on');


%other way around
d = L./D;
%average across recordings
d = squeeze(nanmean(d,1))';

%reorganize by decreasing median
[~,idxsort] = sort(nanmedian(d,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(d(idxsort,:)','Notch','off')
a = get(get(gca,'children'),'children');
t = get(a,'tag'); 
idx = find(strcmp(t,'Box')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', col(i,:),'LineWidth',1.5); 
end
idx = find(strcmp(t,'Median')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', 'k','LineWidth',1.5); 
end
xlabel('Area');
set(gca,'xticklabel',area_label(idxsort),'XTickLabelRotation',45)
ylabel({'Relative dimensionality','(local/subspace)'});
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])

[p,~,stats] = kruskalwallis(d(idxsort,:)',[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
statresults = multcompare(stats,'Display','off');

end


end

