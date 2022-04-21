function Plot_SubspaceOrganization(data)
%Camden MacDowell - timeless
%looks trial-to-trial on the relationship between subspace dimensions
fp = fig_params_cortdynamics;

%get example data from recording 1
[rho,dom,pev,ex] = SubspaceDim_Trial(data{1},0);

[~,area_label] = LoadVariable(data,'rel_performance',[]);

%plot reordered subspace dimension plot
x = sort(ex.rsq,2,'descend');
col = fp.c_area; col = col(strcmp(area_label,'VIS'),:);
figure; hold on;
arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.5,'markersize',4), 1:size(x,1));
plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');  
% shadedErrorBar(1:size(x,2),nanmean(x),std(x,1),'lineprops',{'color',[col,0.25],'linewidth',2});
xlim([0 size(x,2)]); ylim([0 max(x(:))])
xlabel('# of dimensions (sorted)');
ylabel({'Contribution to trial (r^2)'});
title(sprintf('rec 1\n %s | motif %d','VIS',5),'fontweight','normal')
fp.FormatAxes(gca);  box on;
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])

%plot histogram of rhos for the above
figure; hold on; 
histogram(ex.rho(:),'BinWidth',0.2,'FaceAlpha',0.5,'FaceColor',col,'EdgeAlpha',1,'EdgeColor',col)
yval = get(gca,'ylim'); 
plot([nanmean(ex.rho(:)),nanmean(ex.rho(:))],yval,'linewidth',2,'color','k','linestyle','-')
xlabel('rho');
ylabel('# of trial pairs');
fp.FormatAxes(gca); 
box on;
[~,p] = ttest(fisherZ(ex.rho(~isnan(ex.rho))));
title({'Dimension ordering','can differ across trials',sprintf('r=%0.2f p=%0.2f',nanmean(ex.rho(:)),p)},'FontWeight','normal')
fp.FigureSizing(gcf,[5 2 4 4],[2 10 12 12])

%plot average rho across all motifs for example recording (split by area)
col = fp.c_area; 
figure; hold on; 
[~,idxsort] = sort(nanmedian(rho,2),'descend');
boxplot(rho(idxsort,:)','Notch','off')
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
ylabel('Rho');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 14])

[p,~,stats] = kruskalwallis(rho(idxsort,:)',[],'off');
title(sprintf('Rec 1 Areas differ \n in consistency \np=%0.3d',p),'fontweight','normal')
statresults = multcompare(stats,'Display','off');

%example showing the number of times each dimension dominates
figure; hold on;
t=tiledlayout(2,4,'Padding','normal','TileSpacing','normal');
for i = 1:8
    ndim = 5;
    nexttile; hold on; 
    x = squeeze(dom(i,:,1:ndim));
    x = 100*(x./nansum(x,2));
    y = squeeze(pev(i,:,1:ndim));
    y(y<0)=0;
    col = fp.c_area; col = col(strcmp(area_label,area_label{i}),:);
    shadedErrorBar(1:size(x,2),nanmean(x),sem(x,1),'lineprops',{'color',[col,1],'linewidth',2});    
    xlim([0 size(x,2)]); ylim([0 60])
    set(gca,'ytick',[0:20:60])
    xlabel('dimension');
    ylabel({'% of trials'});
    yyaxis right
    %plot the variance
    shadedErrorBar(1:size(y,2),nanmean(y),sem(y,1),'lineprops',{'color',[0.15 0.15 0.15,1],'linewidth',2});   
    ylabel({'Trial r^2'});
    set(gca,'YColor','k');
    title(sprintf('%s',area_label{i}),'fontweight','normal')
    fp.FormatAxes(gca);  box on;    
end
title(t,sprintf('rec 1 | %% of trials dominated by each dimension'),'fontweight','normal')
fp.FigureSizing(gcf,[3 2 5 10],[10 10 20 10])

%% load all data | plot the average Rho per brain area * recording

[rho_all,dom_all,pev_all] = cellfun(@(x) SubspaceDim_Trial(x,0), data,'UniformOutput',0);
%correct for missing loc
temp = rho_all{3};
rho_all{3} = NaN(8,14);
rho_all{3}([1,2,3,5,7,8],:)=temp;

temp = rho_all{4};
rho_all{4} = NaN(8,14);
rho_all{4}([1,2,3,4,5,7,8],:)=temp;


%% compare brain areas
r = cat(3,rho_all{:});
r = reshape(r,size(r,1),size(r,2)*size(r,3));
r = fisherZ(r);

%plot average rho across all motifs split by area)
col = fp.c_area; 
figure; hold on; 
[~,idxsort] = sort(nanmedian(r,2),'descend');
boxplot(r(idxsort,:)','Notch','on')
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
ylabel('Rho_z');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 14])

[p,~,stats] = kruskalwallis(r(idxsort,:)',[],'off');
title(sprintf('Areas differ \n in consistency \np=%0.3d',p),'fontweight','normal')
statresults = multcompare(stats,'Display','off');

%Compare motifs (combine across areas)
figure; hold on; 
r = cat(3,rho_all{:});
r = permute(r,[2,1,3]);
r = reshape(r,size(r,1),size(r,2)*size(r,3));
r = fisherZ(r);
[~,idxsort] = sort(nanmedian(r,2),'descend');
boxplot(r(idxsort,:)','Notch','on','color','k')
a = get(get(gca,'children'),'children');
t = get(a,'tag'); 
idx = find(strcmp(t,'Median')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', 'k','LineWidth',1.5); 
end
xlabel('Motif');
set(gca,'xtick',[1:14],'xticklabel',idxsort,'XTickLabelRotation',90)
ylabel('Rho_z');
fp.FormatAxes(gca);  box on; grid on    
fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 14])
[p,~,stats] = kruskalwallis(r(idxsort,:)',[],'off');
title(sprintf('Motifs Differ in their consistency\n%s p=%0.3f',area_label{cur_a},p),'fontweight','normal')


end %function



% % 
% % 
% % %compare motifs (per brain area) .. not enough samples to actually tests
% % figure; hold on; 
% % for cur_a = 1:8  
% %     subplot(2,4,cur_a); hold on; 
% %     r = cat(3,rho_all{:});
% %     r = squeeze(r(cur_a,:,:));
% %     r = fisherZ(r);
% %     %plot average rho across all motifs split by area)
% %     col = fp.c_area; col = col(strcmp(area_label,area_label{cur_a}),:); 
% %     [~,idxsort] = sort(nanmedian(r,2),'descend');
% %     boxplot(r(idxsort,:)','plotstyle','compact','Colors',col)
% %     xlabel('Motif');
% %     set(gca,'xtick',[1:14],'xticklabel',idxsort,'XTickLabelRotation',90)
% %     ylabel('Rho_z');
% %     fp.FormatAxes(gca);  box on; grid on    
% % 
% %     [p,~,stats] = kruskalwallis(r(idxsort,:)',[],'off');
% %     title(sprintf('%s p=%0.3f',area_label{cur_a},p),'fontweight','normal')
% % end
% % sgtitle(sprintf('Motifs Differ in their consistency'),'fontweight','normal')
% % set(gcf,'units','centimeters','position',[3 3 36 12])













