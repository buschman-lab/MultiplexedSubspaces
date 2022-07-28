function PlotNetworkTSNE(data)
rng('default');
xall = [];
area=[];
for cur_rec = 1:6
%clustering using overlap is included at bottom of page
    x = arrayfun(@(n) reshape(data{cur_rec}(n).rho_all,68*68,10),1:size(data{cur_rec},2),'UniformOutput',0);
    xall{cur_rec} = cat(2,x{:});
    temp = arrayfun(@(n) repmat(data{cur_rec}(n).cur_a,1,10),1:size(data{cur_rec},2),'UniformOutput',0);
    area{cur_rec} = cat(2,temp{:});
    
end
area = cat(2,area{:});
x = cat(2,xall{:});
badidx = sum(isnan(x),2)>1;
x(badidx,:)=[];
% distfun = @(a,b) 1 - ((2*sum((a(:)+b(:))==2))/(sum(a(:)==1)+sum(b(:)==1)));
tcorr_mat = (corr(x));
[cluster_idx, ~, ~] = PhenoCluster(tcorr_mat,'k',12,'louvain_restarts',3,'verbose',1);           
rng('default');
y = tsne(x','distance','correlation','Exaggeration',1,'Perplexity',25);

%% First plot the clustering on TNSE
fp = fig_params_cortdynamics;
figure; hold on; 
idx = 2;
n =numel(unique(cluster_idx(:,idx)));
col = distinguishable_colors(n);
gscatter(y(:,1),y(:,2),cluster_idx(:,idx),col,repmat('.',1,n),ones(1,n)*5)
legend off
fp = fig_params_cortdynamics;
fp.FormatAxes(gca); box on; grid off; set(gca,'xtick',[],'ytick',[]);
fp.FigureSizing(gcf,[3 2 6 6],[2 10 10 10])

% mark our special indices
m = arrayfun(@(n) repmat(data{1}(n).cur_motif,1,10),1:size(data{1},2),'UniformOutput',0);
a = arrayfun(@(n) repmat(data{1}(n).cur_a,1,10),1:size(data{1},2),'UniformOutput',0);
m = cat(2,m{:});
a = cat(2,a{:});
d = repmat(1:10,1,size(data{1},2));
exampleidx = [find((m==8)+(a==2)+(d==6)==3),  find((m==8)+(a==2)+(d==9)==3),...% find((m==11)+(a==7)+(d==5)==3),...
    find((m==6)+(a==4)+(d==5)==3),find((m==6)+(a==4)+(d==10)==3)];

for i = 1:numel(exampleidx)
    plot(y(exampleidx(i),1),y(exampleidx(i),2),'marker','x','markersize',8,'color','k','linestyle','none')
end

exampleidx = [find((m==11)+(a==7)+(d==2)==3),  find((m==6)+(a==4)+(d==3)==3)];
for i = 1:numel(exampleidx)
    plot(y(exampleidx(i),1),y(exampleidx(i),2),'marker','d','markersize',8,'color','k','linestyle','none')
end

%Plot the area colors within each cluster
figure; hold on; 
gscatter(y(:,1),y(:,2),area,fp.c_area,repmat('.',1,8),ones(1,8)*1.5)
legend off
fp = fig_params_cortdynamics;
fp.FormatAxes(gca); box on; grid off; set(gca,'xtick',[],'ytick',[]);
fp.FigureSizing(gcf,[3 2 2.5 2.5],[2 10 10 10])


% compute the number of unique areas per cluster
clust = unique(cluster_idx(:,idx));
n = numel(clust);
grp = NaN(1,n);
for i = 1:n
    grp(i) = numel(unique(area(cluster_idx(:,idx)==clust(i))));
end
figure; hold on; 
histogram(grp,'binwidth',1,'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'facealpha',0.25)
[xboot,stats] = pairedBootstrap(grp',@nanmedian);
title(sprintf('%d clust with 1 region \n %0.0f %0.0f %0.0f',sum(grp==1),nanmedian(xboot),stats.ci(1),stats.ci(2)),'fontweight','normal');
sum(grp==1)
fp.FormatAxes(gca); box on; grid off;
yvals = get(gca,'ylim');
% plot([nanmedian(xboot),nanmedian(xboot)],yvals,'linewidth',1.5,'color','k','linestyle','--')
set(gca,'ylim',yvals);
xlabel('# of regions');
ylabel('# of clusters');
xlim([0.5 8.5])
set(gca,'xtick',[1,4,8])
fp.FigureSizing(gcf,[3 2 2 2.5],[2 10 10 10])

end