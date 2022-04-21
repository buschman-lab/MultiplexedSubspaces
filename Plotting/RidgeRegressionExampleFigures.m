function RidgeRegressionExampleFigures(data,area_name,motif_num)
fp = fig_params_cortdynamics;

%motif_num are the two motifs to compare
[engment,area_label] = LoadVariable(data,'engagement',area_name);
full_mdl = LoadVariable(data,'ridge_performance',area_name);
trial_var = LoadVariable(data,'ttt_variability',area_name);
mean_act = LoadVariable(data,'ttt_activity_mean_all',area_name,[],'none');
mean_act = nanmean(mean_act,3);

col = fp.c_area(ismember(area_label,area_name),:);

m = motif_num(1); 
mout = motif_num(2);

%plot example cross validations Mean +/- SEM for full model 
figure; hold on; 
errorbar(1, nanmean(full_mdl(:,m)),sem(full_mdl(:,m)),'linestyle','none','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
plot(cat(1,ones(1,size(full_mdl,1)),2*ones(1,size(full_mdl,1))),full_mdl(:,[m,mout])','marker','o','markersize',5,'color',[0.65 0.65 0.65],'linestyle','-')
errorbar(2,nanmean(full_mdl(:,mout)),sem(full_mdl(:,mout)),'linestyle','none','marker','o','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
set(gca,'xlim',[0.5 2.5],'xtick',[1,2],'xticklabel',[m,mout])
xlabel('motif'); ylabel('performance (10fold xval)');
title(sprintf('Example Xval Performance motif %d area %s',m,area_name),'fontweight','normal')

%Show how sparse activity is; no baseline subtraction
figure; hold on; 
x = Plot_CompareValueBetweenMotifs(mean_act,[m,mout],fp,'right',col);
ylabel({'Neural Activity','(raw)'}); 
title(sprintf('Max trial averaged\nactivity of %s area',area_name),'fontweight','normal')
ylim([min(x(:))-.1 max(x(:))+0.1])

%'vis' motifs should have more activity in it than, say 'non-visual' motif X
figure; hold on; 
x = Plot_CompareValueBetweenMotifs(engment,[m,mout],fp,'right',col);
ylabel({'Neural Activity','(normalized FR)'}); 
title(sprintf('Trial averaged\nactivity of %s area',area_name),'fontweight','normal')
ylim([min(x(:))-.1 max(x(:))+0.1])

%'vis' motifs should have better performance than, say 'non-visual' motif X
figure; hold on; 
x = Plot_CompareValueBetweenMotifs(full_mdl,[m,mout],fp,'right',col);
ylabel('Performance (r^2)');
title(sprintf('Predicting %s trial-to-trial \n variability from other areas',area_name),'fontweight','normal')
ylim([min(x(:))-.025 max(x(:))+0.025])

%Is this generally true across all recordings and motifs
Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'average',fp,'both','xlabel','Neural Activity','ylabel','Performance','color_flag',0)
title({['Motifs that engage ',area_name],'have better predictive power'},'fontweight','normal');
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])

%Also scales with variance | higher perforamnce
Plot_CorrelateValuesBetweenRecordings(trial_var,full_mdl,'average',fp,'both','xlabel','trial-to-trial \sigma^2','ylabel','Performance','color_flag',0)
title({['Motifs that engage ',area_name],'predict more variance','of more variance'},'fontweight','normal');
fp.FigureSizing(gcf,[3 2 4 4],[10 10 20 10])
end