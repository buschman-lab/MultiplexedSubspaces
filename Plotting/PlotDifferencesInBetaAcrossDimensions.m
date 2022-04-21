function PlotDifferencesInBetaAcrossDimensions(data,area_name, motifs, cur_rec)
if nargin <4; cur_rec = 1; end
fp = fig_params_cortdynamics;
%plot the percent sharing example between motifs across subspace dimensions
%across recordings
[~,area_label] = LoadVariable(data,'rrr_dim',[]);
ndim = 10;

%% Plot the similarity in predictive weightings across dimensions 
%load synaptic weights across dimensions
weights = cell(1,ndim);
for i = 1:ndim
   [weights{i}, area_label] = LoadVariable(data,'rrr_synapticweight',[],i);
end

%get values for your two motifs
y = cellfun(@(x) squeeze(x(cur_rec,motifs(1),strcmp(area_label,area_name),:)),weights,'UniformOutput',0);
y = cat(2,y{:});
x = cellfun(@(x) squeeze(x(cur_rec,motifs(2),strcmp(area_label,area_name),:)),weights,'UniformOutput',0);
x = cat(2,x{:});

%get the difference
xy = abs(y-x);
%split into areas
reg = LoadVariable(data,'beta_region',area_name);
reg = squeeze(reg(cur_rec,motifs(1),:));
xy = arrayfun(@(n) xy(reg==n,:),unique(reg(~isnan(reg))),'UniformOutput',0);
col = fp.c_area; col = col(strcmp(area_label,area_name)==0,:);

figure; hold on; 
for i = 1:size(col,1)
    shadedErrorBar(1:ndim,nanmean(xy{i}),sem(xy{i},1),'lineprops',{'color',[col(i,:),0.25],'linewidth',2});
end
legend(area_label(strcmp(area_label,area_name)==0),'location','bestoutside')
xlim([1 ndim]);
xlabel('# of dimensions');
ylabel('\Delta weighted \beta');
title(sprintf('Betas of \n subspaces for motif %d and %d',motifs(1),motifs(2)),'Fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 3 4],[10 10 15 10])


%See if this approximately matches the magnitude of difference in
%trial-to-trial activity between the two areas. 
ttv = LoadVariable(data,'ttt_activity_mean',[]);
y = squeeze(ttv(cur_rec,motifs(1),strcmp(area_label,area_name)==1,:));
x = squeeze(ttv(cur_rec,motifs(2),strcmp(area_label,area_name)==1,:));
xy = abs(y-x); 
xy = arrayfun(@(n) xy(reg==n),unique(reg(~isnan(reg))),'UniformOutput',0);
temp = cellfun(@(x) nanmean(x), xy);
[~,idx] = sort(temp,'descend');
xy = xy(idx);
figure; hold on; 
for i = 1:size(col,1)    
    bar(i,nanmean(xy{i}),'facecolor',[0.5 0.5 0.5],'EdgeAlpha',0);
    errorbar(i,nanmean(xy{i}),sem(xy{i}),'linewidth',1.5,'Color','k')
end
temp = area_label(strcmp(area_label,area_name)==0);
set(gca,'xtick',1:7,'XTickLabel',temp)



end