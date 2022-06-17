function statresults = PlotSubspaceSharing(data,beta,area_name, motifs,cur_rec,dataout,rsqxval)
if nargin <5; cur_rec = 1; end %example to use
%Camden MacDowell - timeless
fp = fig_params_cortdynamics;
%plot the percent sharing example between motifs across subspace dimensions
%across recordings
[~,area_label] = LoadVariable(data,'rrr_dim',[]);
ndim = 5;
nrec = numel(beta);
cur_a = find(strcmp(area_label,area_name));
y = NaN(nrec,ndim);
for cur_d = 1:ndim
   xbeta= arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
   xbeta = cat(3,xbeta{:});
   y(:,cur_d) = xbeta(motifs(1),motifs(2),:);
end
x=y;
x(x<0)=0;
figure; hold on;
%plot the cross validated performance in the background in blue
xx = cellfun(@(x) squeeze(x(cur_a,motifs(1),:))', rsqxval,'UniformOutput',0);
xx = cat(1,xx{:});
xx(xx<0)=0;
% arrayfun(@(n) plot(1:size(xx,2), xx(n,:),'linestyle',':','marker','none','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',4), 1:size(xx,1));
plot(1:size(xx,2), nanmean(xx),'linestyle',':','color',[0.2 0.2 0.8],'linewidth',1.5,'marker','none');  
arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(x,1));
%plot the example in red
plot(1:size(x,2), x(cur_rec,:),'linestyle','-','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',4,'linestyle','--');
plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);  
xlim([0 size(x,2)]); ylim([0 1]);
xlabel('# of dimensions');
ylabel({'Subspace','Generalization (r^2)'});
title(sprintf('Comparing similarity of \n subspaces for motif %d and %d',motifs(1),motifs(2)),'Fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])


%show that for no motif and no areas does generalization go above X dimensions
temp = arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,:)), 1:size(beta,2),'UniformOutput',0);
temp = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)),temp,'UniformOutput',0);
temp = cat(3,temp{:}); 
temp(temp<0) = NaN; 
temp = sum(isnan(temp),3);
temp = arrayfun(@(n) find(temp(n,:)==6,1,'first')-1, 1:size(temp,1),'UniformOutput',1);
figure; hold on;
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
histogram(temp(:),'BinWidth',1,'EdgeAlpha',0,'FaceColor',col)
xlim([0 10])
fp.FigureSizing(gcf,[3 2 2 2],[10 10 10 10])
ylabel('# of pairs'); xlabel({'# of dimensions','with any generalization'});
title({'Max # of dimensions','with ANY sharing'},'fontweight','normal');
fp.FormatAxes(gca);  box on; grid on

%Total amount of shared performance across dimensions
rrr_mdl = LoadVariable(data,'rel_performance',area_name);

% Get a summary measure of shared dimensioanlity across all dimensions
z = squeeze(rrr_mdl(:,motifs(1),1:ndim));
%convert into rel percent of each dimension
z = [z(:,1),diff(z,[],2)];
x = y; 
x(x<0)=0;
x = cumsum(x.*z,2)*100;
figure; hold on;
arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(x,1));
plot(1:size(x,2), x(cur_rec,:),'linestyle','-','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',4,'linestyle','--');
plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);  
xlim([0 size(x,2)]); ylim([25 75]);
xlabel('# of dimensions');
ylabel({'% of Subspace','that is shared'});
title(sprintf('Comparing similarity of \n subspaces for motif %d and %d',motifs(1),motifs(2)),'Fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])

%% average per recording and organize with dendrogram and visualize for each area[~,area_label] = LoadVariable(data,'rrr_dim',[]);
gen = NaN(numel(area_label),14,14,nrec);
for cur_a = 1:numel(area_label) 
    for m = 1:14
        for mm = 1:14
            y = NaN(nrec,ndim);
            for cur_d = 1:ndim
               xbeta= arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
               xbeta = cat(3,xbeta{:});
               y(:,cur_d) = xbeta(m,mm,:);
            end
            z = squeeze(rrr_mdl(:,m,1:ndim));
            z = ([z(:,1),diff(z,[],2)]); %%%%%%%./z(:,2);%uncomment to convert into rel percent of each dimension, out of the totalfrom all dimensions used 
            x = y; 
            x(x<0)=0;
            gen(cur_a,m,mm,:) = max(cumsum(x.*z,2)*100,[],2);
        end
    end
end


for i = 1:8
    genavg = squeeze(nanmean(gen,4));
    genavg = squeeze(genavg(i,:,:));
    genavg(isnan(genavg))=0;
    Z = linkage(1-genavg,'average','euclidean');
    [~,~,idx] = dendrogram(Z); close; 
    a = genavg(idx,idx);
    a(eye(size(a))==1) = nan;
    figure; hold on; 
    imagesc(a); colorbar
    set(gca,'xtick',1:14,'xticklabel',idx,'ytick',1:14,'yticklabel',idx)
    xlim([0.5 14.5]); ylim([0.5 14.5]); 
    xlabel('Motif');
    ylabel('Motif');
    title(sprintf('%s % Generalization',area_label{i}),'Fontweight','normal')
    fp.FormatAxes(gca);  box on; grid on
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    colormap magma
end


%% same thing, but motifs orgnized by their similarity in pattern of engagement
ttv = LoadVariable(data,'ttt_activity_mean_all',[]);
temp = squeeze(nanmean(ttv,4)); %average within area per recording
temp = squeeze(nanmean(temp,1)); %average across recordings
Z = linkage(temp,'average','euclidean');
% idx = optimalleaforder(Z,pdist(temp));
[~,~,idx] = dendrogram(Z); close; 

for i = 1:8
    genavg = squeeze(nanmean(gen,4));
    genavg = squeeze(genavg(i,:,:));
    genavg(isnan(genavg))=0;
    a = genavg(idx,idx);
    a(eye(size(a))==1) = nan;
    figure; hold on; 
    imagesc(a); colorbar
    set(gca,'xtick',1:14,'xticklabel',idx,'ytick',1:14,'yticklabel',idx)
    xlim([0.5 14.5]); ylim([0.5 14.5]); 
    xlabel('Motif');
    ylabel('Motif');
    title(sprintf('%s % Generalization',area_label{i}),'Fontweight','normal')
    fp.FormatAxes(gca);  box on; grid off
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    colormap magma
end

%% Compare across areas (average per recording so will have # motifs)
genavg = squeeze(nanmean(gen,4));
genavg = reshape(genavg,size(genavg,1),size(genavg,2)*size(genavg,3));
%reorganize by decreasing median
[~,idxsort] = sort(nanmedian(genavg,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(genavg(idxsort,:)','Notch','on')
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
ylabel('% Generalization');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])

[p,~,stats] = kruskalwallis(genavg(idxsort,:)',[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
figure; statresults = multcompare(stats,'Display','on');

%%
%bootstrap distribution
genavg = reshape(gen,size(gen,1),size(gen,2)*size(gen,3)*size(gen,4));
%reorganize by decreasing median
[~,idxsort] = sort(nanmean(genavg,2),'descend');
genboot = pairedBootstrap(genavg',@nanmean);
[~, area_name] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
figure; hold on;       
col = fp.c_area(ismember(area_name,area_label),:);
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
%sort
col = col(idxsort);
CompareViolins(genboot(:,idxsort)',fp,'plotspread',0,'divfactor',0.2,'plotaverage',1,'col',col,'distWidth',0.75);
% if ~isempty(xperm)
%    CompareViolins(xperm(:,ind)',fp,'plotspread',0,'divfactor',spreadfactor(1),'plotaverage',0,'col',repmat({[0 0 0]},1,numel(ind)),'distWidth',spreadfactor(2)); 
% end
set(gca,'XTickLabel',area_label(idxsort),'XTickLabelRotation',45)
yval = get(gca,'ylim');
set(gca,'ylim',[55,yval(2)]);
ylabel('% Generalization');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 3],[10 10 10 10])




%% Do per recording
genavg = squeeze(nanmean(gen,2));
genavg = squeeze(nanmean(genavg,2));
%reorganize by decreasing median
[~,idxsort] = sort(nanmean(genavg,2),'descend');
col = fp.c_area; 
figure; hold on; 
plot(nanmean(genavg(idxsort,:)'),'linewidth',1.5,'color',[0.75 0.75 0.75]);
arrayfun(@(n) errorbar((n),nanmean(genavg(idxsort(n),:)'),sem(genavg(idxsort(n),:)'),'color',col(n,:),'linewidth',2,'marker','o','markersize',fp.markersizesmall),1:numel(idxsort));
% arrayfun(@(n) plot(genavg(idxsort,n),'color',[0.75 0.75 0.75],'linewidth',1),1:size(genavg,2));
xlabel('Area');
set(gca,'xtick',1:8,'xticklabel',area_label(idxsort),'XTickLabelRotation',45)
ylabel('% Generalization');
title({'Averaged across motifs','per recording'},'fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 10 10])

%% compare with the number of dimensions that region has in it's subspace?
rrr_d = LoadVariable(data,'rrr_dim',[]);
rrr_d = squeeze(nanmean(rrr_d,1))';
%reorganize by decreasing median
[~,idxsort] = sort(nanmedian(rrr_d,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(rrr_d(idxsort,:)','Notch','off')
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
ylabel('Subspace Dimensionality');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 10 12])
title({'% Generalization is not','due to dimensionality of subspace'},'fontweight','normal');

%% What is driving this? One hypothesis is that it's because a region is more of an integrator or segregator
%compare V in with B out

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


% Correlated with your generalization index
genavg = squeeze(nanmean(gen,4));
genavg = reshape(genavg,size(genavg,1),size(genavg,2)*size(genavg,3));
%reorganize by decreasing median
a = nanmedian(genavg,2);
b = nanmedian(x,2);
col = fp.c_area; 
Plot_CorrelateValuesBetweenRecordings(a,b,'combo',fp,'right','xlabel',{'% Generalization'},...
    'ylabel',{'In-Out Representational','Similarity'},'color_flag',1,'corrtype','spearman','col',col,'addjitter',0)
xlim([min(a(:))-2.5 max(a(:))+2.5])
title(sprintf('Similarity Underlies Generalizability'),'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%% Now do in the opposite direction
y = NaN(nrec,cur_a,14,ndim);
for cur_a = 1:8
    for cur_d = 1:ndim
        LocalIn = LoadVariable(data,'rrr_beta',area_label(cur_a),cur_d);
        LocalOut = LoadVariable(dataout,'rrr_V',area_label(cur_a),cur_d);
        for cur_rec = 1:6
            for cur_m = 1:14
                a = squeeze(LocalIn(cur_rec,cur_m,:));
                b = squeeze(LocalOut(cur_rec,cur_m,:));
                y(cur_rec,cur_a,cur_m,cur_d) = corr(a,b,'rows','complete');
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

% Correlated with your generalization index
genavg = squeeze(nanmean(gen,4));
genavg = reshape(genavg,size(genavg,1),size(genavg,2)*size(genavg,3));
%reorganize by decreasing median
a = nanmedian(genavg,2);
b = nanmedian(x,2);
col = fp.c_area; 
Plot_CorrelateValuesBetweenRecordings(a,b,'combo',fp,'right','xlabel',{'% Generalization'},...
    'ylabel',{'In-Out Representational','Similarity'},'color_flag',1,'corrtype','spearman','col',col,'addjitter',0)
xlim([min(a(:))-2.5 max(a(:))+2.5])
title(sprintf('Similarity Underlies Generalizability'),'fontweight','normal') 
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])



%% Plot showing that effect was dominated at lower dimensions
y = NaN(nrec,numel(area_label),14,14,ndim);
for cur_a = 1:numel(area_label)
    for m = 1:14
        for mm = 1:14
            for cur_d = 1:ndim
                xbeta= arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
                xbeta = cat(3,xbeta{:});
                y(:,cur_a,m,mm,cur_d) = xbeta(m,mm,:);
            end
        end
    end
end
y = squeeze(nanmean(y,1))*100;
figure; hold on;
for cur_a = 1:8
    temp = squeeze(y(cur_a,:,:,:));
    temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
    temp(temp<0)=0;
    shadedErrorBar(1:ndim,nanmean(temp),sem(temp,1),'lineprops',{'color',[col(cur_a,:),0.25],'linewidth',2});
end
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 4 4],[10 10 10 12])
ylabel('% Generalization');
title({'% Generalization is','due to lower dimensions'},'fontweight','normal');
xlabel('Subspace dimensions')




end %function end
