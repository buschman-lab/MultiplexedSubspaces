function PlotSubspaceBidirectionality(data_nolag,data,dataout_nolag,dataout)
%Camden MacDowell - timeless
fp = fig_params_cortdynamics;

%Compare the total explained variance
in = squeeze(nanmean(LoadVariable(data,'ridge_performance',[]),1));
out = squeeze(nanmean(LoadVariable(dataout,'ridge_performance',[]),1));
innolag = squeeze(nanmean(LoadVariable(data_nolag,'ridge_performance',[]),1));
outnolag = squeeze(nanmean(LoadVariable(dataout_nolag,'ridge_performance',[]),1));

%plot the correlation between the two across all motifs and areas
col = fp.c_area; 
Plot_CorrelateValuesBetweenRecordings(in',out','combo',fp,'right','xlabel',{'% Generalization'},...
    'ylabel',{'In-Out Representational','Similarity'},'color_flag',1,'corrtype','spearman','col',col,'addjitter',0)

%are any of them more predictive with a lag than without
y = in-innolag; 
x = out-outnolag; 



























end 



ndim = 10;
y = NaN(nrec,cur_a,14,ndim);
for cur_a = 1:8
    for cur_d = 1:ndim
        LocalIn = LoadVariable(data,'rrr_V',area_label(cur_a),cur_d);
        LocalOut = LoadVariable(dataout,'rrr_beta',area_label(cur_a),cur_d);
        for cur_rec = 1:6
            for cur_m = 1:14
                a = squeeze(LocalIn(cur_rec,cur_m,:));
                b = squeeze(LocalOut(cur_rec,cur_m,:));
                y(cur_rec,cur_a,cur_m,cur_d) = corr(a,b,'rows','complete','type','pearson'); 
            end
        end
    end
end
x = squeeze(nanmean(nanmean(fisherZ(y),1),4));
[~,idxsort] = sort(nanmedian(x,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(x(idxsort,:)','Notch','on')
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
ylabel('Integration');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 10 12])
[p,~,stats] = kruskalwallis(x(idxsort,:)',[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
figure; statresults = multcompare(stats,'Display','on');