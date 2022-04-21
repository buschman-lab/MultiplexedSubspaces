function AnalyzeCCs_draft4()
% camden - timeless

%Load the data_rrr for our specific motifs from each rec
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';

data = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDm\w*.mat'],0,{folder}); 
%     [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPm\w*.mat'],0,{folder}); 
    data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','inactive_idx','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif'),fn);
end

datacca = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*CCA_muaflag1_GROUPED\w*.mat'],0,{folder}); 
    datacca{cur_rec} = cellfun(@(x) load(x,'r','area_label','motif'),fn);
end

datar = cell(1,6); %reversed
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDREVERSE\w*.mat'],0,{folder}); 
    datar{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif'),fn);
end


%remove noise motif and null motif
for i = 1:6
   idx = ismember(arrayfun(@(n) datacca{i}(n).motif, 1:size(datacca{i},2)),[2,16]);
   datacca{i}(idx)=[];     
   idx = ismember(arrayfun(@(n) datar{i}(n).motif, 1:size(datar{i},2)),[2,16]);
   datar{i}(idx)=[];      
   idx = ismember(arrayfun(@(n) data{i}(n).motif, 1:size(data{i},2)),[2,16]);
   data{i}(idx)=[];     
end


%% Question 1: When brain area is engaged by a motif, does it have more predictable varaibilty with other areas in the brian? 
% Look at visual area
fp = fig_params_cortdynamics;
area_name = 'VIS';
engment = LoadVariable(data,'engagement',area_name);
full_mdl = LoadVariable(data,'ridge_performance',area_name);
trial_var = LoadVariable(data,'ttt_variability',area_name);
m = 5;
mout = 6;

%plot example recording Mean +/- SEM for full model 
figure; hold on; 
errorbar(1, nanmean(full_mdl(:,m)),sem(full_mdl(:,m)),'linestyle','none','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
plot(cat(1,ones(1,size(full_mdl,1)),2*ones(1,size(full_mdl,1))),full_mdl(:,[m,mout])','marker','o','markersize',5,'color',[0.65 0.65 0.65],'linestyle','-')
errorbar(2,nanmean(full_mdl(:,mout)),sem(full_mdl(:,mout)),'linestyle','none','marker','o','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
set(gca,'xlim',[0.5 2.5],'ylim',[0.07 0.13],'xtick',[1,2],'xticklabel',[m,mout])
xlabel('motif'); ylabel('performance (10fold xval)');
title(sprintf('Example Performance motif %d area %s',m,area_name),'fontweight','normal')

%'vis' motifs should have more activity in it than, say 'non-visual' motif X
figure; hold on; 
Plot_CompareValueBetweenMotifs(engment,[m,mout],fp,'right');
ylabel({'Neural Activity','(normalized FR)'}); 
title(sprintf('Trial averaged\nactivity of %s area',area_name),'fontweight','normal')

%'vis' motifs should have better performance than, say 'non-visual' motif X
figure; hold on; 
Plot_CompareValueBetweenMotifs(full_mdl,[m,mout],fp,'right');
ylabel('Performance (r^2)');
title(sprintf('Predicting %s trial-to-trial \n variability from other areas',area_name),'fontweight','normal')
ylim([0.06 0.16])

%Is this generally true across all recordings and motifs
Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'zscore',fp,'both','xlabel','Neural Activity','ylabel','Performance','color_flag',0)
title({'zscore within each rec'});

%controls showing that motif 5 is high up on the positive side of this (as expected)
% Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'average',fp,'both','xlabel','Neural Activity','ylabel','Performance','color_flag',0,'focus_motifs',[m,mout])
% title({'average across motifs'});
% Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'maxdiv',fp,'both','xlabel','Neural Activity','ylabel','Performance','color_flag',0,'focus_motifs',m)
% title({'Normalized by Max'})
% Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'none',fp,'both','xlabel','Neural Activity','ylabel','Performance','color_flag',2)
% title({'no normalizeation normalized'})
% 
% %control analysis: this is why we cannot directly combine across recordings
% Plot_CorrelateValuesBetweenRecordings(engment,full_mdl,'combo',fp,'both','xlabel','Neural Activity','ylabel','Performance')
% title({'performance and activity','systeamically vary across recs'},'fontweight','normal')

%also it scales with the variance, so it's explaining a higher amount of
%stuff with higher variance (so we should weight by that later)
Plot_CorrelateValuesBetweenRecordings(trial_var,full_mdl,'average',fp,'both','xlabel','trial-to-trial \sigma^2','ylabel','Performance','color_flag',0)
title({['Motifs that engage ',area_name],'predict more variance','of more variance'},'fontweight','normal');

%so we should scale performance by that value
figure; hold on; 
Plot_CompareValueBetweenMotifs(full_mdl.*trial_var,[m,mout],fp,'right');
ylabel({'Weighted','performance'});
title(sprintf('Predicting %s trial-to-trial \n variability from other areas (scaled)',area_name),'fontweight','normal')

%similarly, you can plot an example of two regions engaged with visual cortex
m = [5,7];
mout =6;

%'vis' motifs should have better performance than, say 'non-visual' motif X
figure; hold on; 
Plot_CompareValueBetweenMotifs(full_mdl.*trial_var,[m,mout],fp,'right');
ylabel({'Weighted','performance'});
title(sprintf('Predicting %s trial-to-trial \n variability from other areas',area_name),'fontweight','normal')
ylim([1 7])
% handles = findall(groot,'Type','figure');
% saveCurFigs(handles,{'-dpng','-dsvg'},'FullModel',savedir,0); close all    

%% Test the hypothesis that interactions confined to subspaces
%plot example recording 
m = 5;
area_name = 'VIS';
rrr_mdl = LoadVariable(data,'rel_performance',area_name);
ndim = 15;
thresh = 0.75;
figure; hold on;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
rrr_mdl = squeeze(rrr_mdl(:,m,1:30));
arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(rrr_mdl,1));
plot(1:size(rrr_mdl,2), nanmean(rrr_mdl),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);  
d = arrayfun(@(n) find(rrr_mdl(n,:)>=thresh,1,'first'),1:size(rrr_mdl,1));
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel({'Performance','(r^2 of explainable variance)'});
title(sprintf('Motif %d needed %d to %d dim \n to capture %0.1d%% PEV',m,min(d), max(d),thresh*100),'Fontweight','normal')
fp.FormatAxes(gca);  box on;
fp.FigureSizing(gcf,[3 2 6 6],[10 10 10 10])

%This was not due to our regression approach.
figure; hold on; 
cca_d = LoadVariable(datacca,'cca_dim',area_name);
rrr_d = LoadVariable(data,'rrr_dim',area_name);
Plot_CompareValueBetweenMotifs([rrr_d(:,m),cca_d(:,m)],[1,2],fp,'both');
ylabel('# of dimensions'); xlabel('method'); set(gca,'xtick',[1,2],'xticklabel',{'RRR','CCA'})
title(sprintf('Low dimensionality \nis not due to anlaysis \nmethod | Motif %d',m),'Fontweight','normal')
ylim([0 12])


%Relationship between local dimensionality and the subspace dimensionality
rrr_d = LoadVariable(data,'rrr_dim',area_name);
% localdim = LoadVariable(data,'svca',area_name,0.05);
localdim = LoadVariable(data,'pca_auc_time',area_name,80);
Plot_CorrelateValuesBetweenRecordings(rrr_d,localdim,'combo',fp,'right','xlabel',{'Subspace Dimensionality'},'ylabel','Local dimensionality','color_flag',2,'corrtype','pearson')
plot([0 max(cat(1,rrr_d(:),localdim(:)))],[0 max(cat(1,rrr_d(:),localdim(:)))],'linestyle',':','color','k')
title({'Subspace dimensionality is','lower than local dimensionality'},'fontweight','normal') 


%Nor how we measured dimensionality | correlation between dimensionality cut off and AUC
area_name = 'VIS';
d = LoadVariable(data,'rrr_dim',area_name);
auc = LoadVariable(data,'rrr_auc',area_name,15);
Plot_CorrelateValuesBetweenRecordings(d,auc,'average',fp,'both','xlabel',{'# of dimenions to explain','80% of variance'},'ylabel','1-AUC','color_flag',0,'corrtype','kendall')
title({'Using AUC as a','measure of dimensionality'},'fontweight','normal') 

handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng','-dsvg'},'CommSubspace',savedir,0); close all    

%same test for local
% localdim = LoadVariable(data,'svca',area_name,0.05);
% Plot_CorrelateValuesBetweenRecordings(localdim,localauc,'combo',fp,'both','xlabel',{'# of dimenions to explain','80% of variance'},'ylabel','1-AUC','color_flag',0,'corrtype','spearman')


%% Do these regression models make sense (i.e., synaptic weights follow expected distributions)
%in a motif that engages RSP and VIS is the prediction of VIS weighted by RSP 
area_name = 'VIS';
m = 5;
mout = 6;
b = LoadVariable(data,'rrr_beta',area_name,1);
[reg,area_lbl] = LoadVariable(data,'beta_region',area_name);
fh = PlotRegressionWeights(b,m,reg,area_lbl,fp);
title(sprintf('raw beta weights predicting %s activity during motif %d',area_name,m),'fontweight','normal');
xlabel('neurons'); ylabel('Beta Weights');

%weight by firing rate
b = LoadVariable(data,'rrr_synapticweight',area_name,1);
fh = PlotRegressionWeights(b,m,reg,area_lbl,fp);
title(sprintf('synaptic weights into %s activity during motif %d',area_name,m),'fontweight','normal');
xlabel('neurons'); ylabel('Synaptic Weights');

%Plot ordered per area
b = LoadVariable(data,'rrr_synapticweight_sorted',area_name,1);
fh = PlotRegressionWeights(b,m,reg,area_lbl,fp);
title(sprintf('synaptic weights into %s activity during motif %d',area_name,m),'fontweight','normal');
xlabel('neurons'); ylabel('Synaptic Weights');

%compare the synaptic weights between a motif that engages the region and one that doesn't
[b,area_lbl] = LoadVariable(data,'synapticweight_sum_abs',[],1);
figure; hold on;
PlotRegressionWeightsByArea(b,m,mout,area_lbl,area_name,fp,[2,2,2,2,2,2,2])

handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng','-dsvg'},'Betavalidationm5and6',savedir,0); close all    

% Together this supports the idea of communication subspaces as a
% generalized feature of neural populations dynamics



%% Are these subspaces shared 


























%% are these subspaces the same? or different
[rsq,matchdim] = SharedBetas(data{1},'B',0);
[rsq_shuf,~] = SharedBetas(data{1},'B',1);

%% are these subspaces relying on different neurons? or just increase/decrease weights on individual neurons. 
%first look at dimension 1
% n_auc = LoadVariable(data,'dom_neurons_auc',area_name,1);
trial_var = LoadVariable(data,'ttt_variability',area_name);
for cur_d = [1,2,3,5,10]
    b = LoadVariable(data,'rrr_synapticweight',area_name,cur_d);
    rsq = LoadVariable(data,'rel_performance',area_name);
    rsq = rsq(:,:,cur_d);%.*trial_var;

    %get the correlation between the synaptic weights across motifs
    [rho,~,idx] = CorrelatedNeuronContributions(b,rsq);
    figure; hold on; 
    imagesc(nanmean(rho,3)); colormap magma; c=colorbar;
    xlim([0.5 14.5]); ylim([0.5 14.5]); xlabel('motif'); ylabel('motif');
    set(gca,'xtick',1:14,'xticklabel',idx(1:14))
    set(gca,'ytick',1:14,'yticklabel',idx(1:14))
    ylabel(c,'rho')
    fp.FigureSizing(gcf,[2 2 6 6],[10 10 12 12]);
    fp.FormatAxes(gca);
    title(sprintf('Correlation in Betas Dimension %d',cur_d),'fontweight','normal')


    figure; plot(nanmean(rsq(:,idx)))
    fp.FormatAxes(gca);
    title(sprintf('Strength of subspace %d',cur_d),'fontweight','normal')
    set(gca,'xtick',1:14,'xticklabel',idx(1:14))
    xlabel('motif'); ylabel({'weighted','performance'});
    fp.FigureSizing(gcf,[2 2 6 2],[10 10 12 6]);
end

% %get the correlation between first dimension performance and the number of neurons
% Plot_CorrelateValuesBetweenRecordings(rsq,n_auc,'average',fp,'left','xlabel','weighted performance','ylabel','synaptic dimensionality','color_flag',0,'corrtype','pearson')
% 
% %Control: importantly, all dimensions in all motifs required a fair number of neurons
% n_dom = LoadVariable(data,'dom_neurons_num',area_name,1);
% fract_dom = n_dom./size(b,3)*100;

%in otherwords, are the weaker distributed brain signals that are
%'constantly' there due large and variable networks that are not shared
%across motifs.
%however the relationship goes away at higher dimensions, showing how these
%are indeed different patterns. 

% handles = findall(groot,'Type','figure');
% saveCurFigs(handles,{'-dpng','-dsvg'},'VIS IN subspace betas',savedir,0); close all    

%% Easy, straightforward plot of reciprocity
%need to show that some areas that actively engage a region have
%significant directionality. 
area_name = 'PRE';
pin = LoadVariable(data,'ridge_performance',area_name);
pout = LoadVariable(datar,'ridge_performance',area_name);
ain = LoadVariable(data,'rrr_dim',area_name);
aout = LoadVariable(datar,'rrr_dim',area_name);
pin = 100*pin;
pout = 100*pout;

din = pin./ain;
dout = pout./aout;

d = din+dout;%din+dout;
r = pin+pout; %pin+pout;

%individual in and out is not correlated with performance
Plot_CorrelateValuesBetweenRecordings(pin,ain,'average',fp,'both','xlabel','reciprocity','ylabel','subspace strength','color_flag',0,'corrtype','pearson');

%% how recipricol versus directional a brain 
%plot directionality of motifs
figure; hold on; 
m=[5,7,8]
Plot_CompareValueBetweenMotifs(directional,m,fp,'right');
ylabel({'Weighted','performance'});




%% 
%individual in and out is not correlated with performance
Plot_CorrelateValuesBetweenRecordings(r,d,'average',fp,'both','xlabel','reciprocity','ylabel','subspace strength','color_flag',0,'corrtype','pearson')


%also plot directionality... show is different between diff motifs (as one
%might expect) and different between different areas. 

%correlation between performance in and out
Plot_CorrelateValuesBetweenRecordings(pin,pout,'average',fp,'both','xlabel','performance in','ylabel','performance out','color_flag',0,'corrtype','pearson')

%correlation between dimensionality in and out
Plot_CorrelateValuesBetweenRecordings(ain,aout,'average',fp,'both','xlabel','# dim in','ylabel','# dim in','color_flag',0,'corrtype','kendall')

%correlation between subspace in and out strength (i.e. effective dimensionality)
Plot_CorrelateValuesBetweenRecordings(din,dout,'average',fp,'both','xlabel','subspace strength (in)','ylabel','subspace strength (out)','color_flag',0,'corrtype','pearson')

%individual in and out is not correlated with performance
Plot_CorrelateValuesBetweenRecordings(r,d,'average',fp,'both','xlabel','reciprocity','ylabel','subspace strength','color_flag',0,'corrtype','pearson')

%same for CCA
cca_d = LoadVariable(datacca,'cca_dim',area_name);
cca_r = LoadVariable(datacca,'cca_strength',area_name);
temp = cca_r./cca_d;
Plot_CorrelateValuesBetweenRecordings(r,temp,'average',fp,'both','xlabel','performance','ylabel','bang for buck','color_flag',0,'corrtype','pearson')








% 
% 
% 
% 
% 
% 
% poutc = pout; 
% for i = 1:size(pin,1)
%     coef = polyfit(pout(i,:),pin(i,:),1); 
%     correction = @(x) coef(1)*x+coef(2);
%     poutc(i,:) = correction(pout(i,:));
% end
% r =pin+poutc;
% directional = (pin-poutc)./(pin+poutc); %-1 = out>in. 0=out=in  1= in>out
% 
% 
% doutc = dout; 
% for i = 1:size(pin,1)
%     coef = polyfit(dout(i,:),din(i,:),1); 
%     correction = @(x) coef(1)*x+coef(2);
%     doutc(i,:) = correction(dout(i,:));
% end
% d = doutc+din;
% 
























%% What could be underlying these relationships
[d_time,area_label] = LoadVariable(data,'pca_auc_time',[],80);
[d_time,area_label] = LoadVariable(data,'svca_trialvar_acrosstime',[]);
rin_time = LoadVariable(data,'rsq_acrosstime',[],15);
ain_time = LoadVariable(data,'auc_acrosstime',[],15);
tavg = LoadVariable(data,'trialavg_timecourse',[]);
tvar = LoadVariable(data,'trialvar_timecourse',[]);
tvar_r = LoadVariable(datar,'trialvar_timecourse_rev',[]);
rout_time = LoadVariable(datar,'rsq_acrosstime_rev',[],15);
aout_time = LoadVariable(datar,'auc_acrosstime_rev',[],15);
trace = LoadVariable(data,'trialavg_timecourse',[]);
trace_other = LoadVariable(data,'trialavg_timecourse_other',[]);
trial_var = LoadVariable(data,'ttt_variability',[]);
trial_activity = LoadVariable(data,'engagement',[]);


%% Plot example of neural activity for motif 5 over time
m = 7;
area_name = 'VIS';
figure; hold on;
Plot_VariableOverTime(squeeze(rin_time(:,m,strcmp(area_label,area_name),:)).*squeeze(tvar(:,m,strcmp(area_label,area_name),:)),fp.c_ff,'none');
ylabel('performance')
yyaxis right
Plot_VariableOverTime(squeeze(ain_time(:,m,strcmp(area_label,area_name),:)),fp.c_lr,'minmax');
% Plot_VariableOverTime(squeeze(d_time(:,m,strcmp(area_label,area_name),:)),fp.c_glm,'minmax');
hold on; 
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])
legend({'IN predictive power','Dimensionality'})
title('Predictive power and AUC over time','fontweight','normal')

figure; hold on;
Plot_VariableOverTime(squeeze(rout_time(:,m,strcmp(area_label,area_name),:)).*squeeze(tvar_r(:,m,strcmp(area_label,area_name),:)),fp.c_ff,'minmax');
ylabel('performance')
yyaxis right
Plot_VariableOverTime(squeeze(aout_time(:,m,strcmp(area_label,area_name),:)),fp.c_lr,'minmax');
% Plot_VariableOverTime(squeeze(d_time(:,m,strcmp(area_label,area_name),:)),fp.c_glm,'minmax');
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])
title('Predictive power and AUC over time','fontweight','normal')
legend({'OUT predictive power','Dimensionality'})


handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng','-dsvg'},'Motif 5 VIS Over time',savedir,0); close all    


%% Plot Example 
%plot the effective dimensionality (flipped) versus the trial average
%the plot the
m = 10;
area_name = 'VIS';
figure; hold on;
x = rin_time./ain_time;
y = rout_time./aout_time;
Plot_VariableOverTime(squeeze(x(:,m,strcmp(area_label,area_name),:)),fp.c_ff,'none');
% Plot_VariableOverTime(squeeze(y(:,m,strcmp(area_label,area_name),:)),fp.c_lr,'minmax');
% ylim([-1 1]); 
set(gca,'ydir','reverse')
ylabel('effective dimensionality')
%plot the trial average activity
yyaxis right; 
Plot_VariableOverTime(squeeze(trace_other(:,m,strcmp(area_label,area_name),:)),fp.c_none,'minmax');
ylabel('neural activity')
xlabel('time from motif onset (ms)');
legend({'In subspace','Out Subspace','Trial average'})
set(gca,'YColor','k')
fp.FormatAxes(gca)
xlim([1,10])

% ylim([-1 1]);
title(sprintf('Neural populations settle into lower \n dimensional subspace when interacting'),'fontweight','normal')


%% Look at effective dimensionality and reciprocity in the visual cortex
area_name = 'VIS';
m = [7,5];
x = rin_time.*tvar./ain_time;
x = squeeze(x(:,:,strcmp(area_label,area_name),:));
y = rout_time.*tvar_r./aout_time;
y=squeeze(y(:,:,strcmp(area_label,area_name),:));

[~,~,rho] = PerformanceXCorr(x,y,1,0);
figure; hold on;
rho = fisherZ(rho);
Plot_CompareValueBetweenMotifs(rho,m,fp,'right');
title('reciprocity','fontweight','normal')

figure; hold on;
temp = -1*max(x.*y,[],3);
Plot_CompareValueBetweenMotifs(temp,m,fp,'left');
title('effective dimensionality','fontweight','normal')

handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng','-dsvg'},'VisualMotifsDimensionalityAndReciprocity',savedir,0); close all    



%% Look at effective dimensionality across all areas

%show a paired test for visual motif 5 and 7
close all
%ward is good for VIS and MOs
id = [1,1,1,0,0,0,0,1];
for i = 1:8
area_name = area_label{i};
% figure; hold on;
rng('default')
Z = linkage([nanmean(squeeze(trial_var(:,:,strcmp(area_label,area_name))))',nanmean(squeeze(trial_activity(:,:,strcmp(area_label,area_name))))'],'ward');
% dendrogram(Z)
clustID = cluster(Z,'MaxClust',2)==1;

mstrong = find(clustID==id(i));
% if i==1
%     mstrong = [1,2,5,6,8,9,10,12];
% end
m = 1:14; 
% m=[5,7,8,12];
% x = rin_time.*tvar./ain_time;
x = rin_time./ain_time;
x = squeeze(x(:,:,strcmp(area_label,area_name),:));
% y = rout_time.*tvar_r./aout_time;
y = rout_time./aout_time;
y=squeeze(y(:,:,strcmp(area_label,area_name),:));
[xrho,lag,rho] = PerformanceXCorr(x,y,1,0);
% figure; hold on;
% rho = fisherZ(rho);
% Plot_CompareValueBetweenMotifs(rho,[m],fp,'left');

temp = -1*max(x,[],3).*max(y,[],3);
% figure; hold on;
% temp = max(x.*y,[],3);
% Plot_CompareValueBetweenMotifs(temp,[m],fp,'left');
%plot correlation versus their peak activity
% arrayfun(@(n) corr(rho(n,m)',temp(n,m)'),1:size(temp,1),'UniformOutput',1)
[~,rho,pval,xtemp,ytemp] = Plot_CorrelateValuesBetweenRecordings(rho(:,m),temp(:,m),'zscore',fp,'left','xlabel','Reciprocity','ylabel','Strength','color_flag',0,'corrtype','pearson');
xtemp = xtemp(:,mstrong);
ytemp = ytemp(:,mstrong);
h = plot(fitlm(xtemp(:),ytemp(:))); 
fp.FormatLSline(h);
delete(h(3:4))
xlabel('Reciprocity');
ylabel('Effective Dimensionality');
h(2).Color=fp.c_ff;
set(gca,'xlim',[-3 3]); set(gca,'ylim',[-3 3])
title(sprintf('%s',area_label{i}))
hh = line(nan, nan, 'Color', 'none');
legend(hh, {sprintf('r=%0.3f \n p=%0.3f',rho,pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best') 

end

% % savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Subspaces'
% handles = findall(groot,'Type','figure');
% saveCurFigs(handles,{'-dpng','-dsvg'},'effectiveDimvsreciprocity',savedir,0); close all    












%% DONE













% really what we want is 'bits' per dimension
m = 5;
area_name = 'VIS';
figure; hold on;
% Plot_VariableOverTime((squeeze(rin_time(:,m,:)))./squeeze(ain_time(:,m,:)),fp.c_ff,'none');
% yyaxis right
% Plot_VariableOverTime((squeeze(rout_time(:,m,:)))./squeeze(aout_time(:,m,:)),fp.c_lr,'none');
Plot_VariableOverTime((squeeze(rin_time(:,m,:)).*squeeze(tvar(:,m,:)))./squeeze(ain_time(:,m,:)),fp.c_ff,'minmax');
yyaxis right
Plot_VariableOverTime((squeeze(rout_time(:,m,:)).*squeeze(tvar_r(:,m,:)))./squeeze(aout_time(:,m,:)),fp.c_lr,'minmax');
% Plot_VariableOverTime((squeeze(rin_time(:,m,:)).*squeeze(tvar(:,m,:))),fp.c_ff,'none');
% yyaxis right
% Plot_VariableOverTime((squeeze(rout_time(:,m,:)).*squeeze(tvar_r(:,m,:))),fp.c_lr,'none');
% Plot_VariableOverTime((squeeze(tvar(:,m,:)))./squeeze(d_time(:,m,:)),fp.c_glm,'none');
% yyaxis right
% Plot_VariableOverTime((squeeze(tvar_r(:,m,:)))./squeeze(d_time(:,m,:)),fp.c_none,'none');

%plot the trial average 
% Plot_VariableOverTime(squeeze(trace(:,m,:)),fp.c_none,'minmax');
% Plot_VariableOverTime(squeeze(trace_other(:,m,:)),fp.c_none,'minmax');
ylabel('performance')

%those with more coordinated sharing, also have stronger subspace (given by
%effective dimensionality
m = find(T==1);
% m = [5,7,8,12];
x = rin_time.*tvar./ain_time;
y = rout_time.*tvar_r./aout_time;
[xrho,lag,rho] = PerformanceXCorr(x,y,1,1);
figure; hold on;
rho = fisherZ(rho);
Plot_CompareValueBetweenMotifs(rho,[m],fp,'left');

x = rin_time.*tvar./ain_time;
y = rout_time.*tvar_r./aout_time;
temp = max((x+y)./2,[],3);
figure; hold on;
% temp = max(x.*y,[],3);
Plot_CompareValueBetweenMotifs(temp,[m],fp,'left');
%plot correlation versus their peak activity
% arrayfun(@(n) corr(rho(n,m)',temp(n,m)'),1:size(temp,1),'UniformOutput',1)
Plot_CorrelateValuesBetweenRecordings(rho(:,m),temp(:,m),'average',fp,'left','xlabel','Reciprocity','ylabel','Strength','color_flag',0,'corrtype','pearson')



%% Reciprocity is the covariation
[rho,lag] = PerformanceXCorr(rin_time(:,:,:),rout_time(:,:,:),1,1);
rho = fisherZ(rho);


figure; hold on;
Plot_CompareValueBetweenMotifs(rho,[5,7],fp,'left');


%%
area_name = 'VIS';
ain = LoadVariable(data,'rrr_auc',area_name,15);
pin = LoadVariable(data,'ridge_performance',area_name);
aout = LoadVariable(datar,'rrr_auc',area_name,15);
pout = LoadVariable(datar,'ridge_performance',area_name);

poutc = pout; 
for i = 1:size(pin,1)
    coef = polyfit(pout(i,:),pin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    poutc(i,:) = correction(pout(i,:));
end
r =pin+pout;

aoutc = aout; 
for i = 1:size(ain,1)
    coef = polyfit(aout(i,:),ain(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aoutc(i,:) = correction(aout(i,:));
end
d = ain+aoutc;


%%
area_name = 'VIS';



%get the correlation between subspace and the activity
%cross correlation between in and out performance, per motif
[rho,lag] = PerformanceXCorr(rin_time,rout_time,1,1);

%get the correlation between the subspace and the local dimensionality
[rho,lag] = PerformanceXCorr(ain_time,d_time,1,1);
%compare with motifs directionality
Plot_CorrelateValuesBetweenRecordings(rho,lag,'average',fp,'left','xlabel','weighted performance','ylabel','synaptic dimensionality','color_flag',0,'corrtype','pearson')









%%

area_name = 'VIS';
localdim = LoadVariable(data,'svca',area_name,0.05);
ain = LoadVariable(data,'rrr_auc',area_name,15);
pin = LoadVariable(data,'ridge_performance',area_name);
aout = LoadVariable(datar,'rrr_auc',area_name,15);
pout = LoadVariable(datar,'ridge_performance',area_name);

poutc = pout; 
for i = 1:size(pin,1)
    coef = polyfit(pout(i,:),pin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    poutc(i,:) = correction(pout(i,:));
end
r = (pin+poutc);

aoutc = aout; 
for i = 1:size(ain,1)
    coef = polyfit(aout(i,:),ain(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aoutc(i,:) = correction(aout(i,:));
end
d = (ain+aoutc)%./localdim;


%%
%reciprocity
Plot_CorrelateValuesBetweenRecordings(r,d,'average',fp,'right','xlabel','reciprocity','ylabel','dimensionality','color_flag',0,'corrtype','pearson')

%separately
Plot_CorrelateValuesBetweenRecordings(pin,ain,'average',fp,'right','xlabel','reciprocity','ylabel','dimensionality','color_flag',0,'corrtype','pearson')
Plot_CorrelateValuesBetweenRecordings(pout,aout,'average',fp,'left','xlabel','reciprocity','ylabel','dimensionality','color_flag',0,'corrtype','pearson')

%in versus out
x = abs(pin-poutc)*100;
y = ain+aoutc;
Plot_CorrelateValuesBetweenRecordings(x,y,'average',fp,'right','xlabel','directionality','ylabel','dimensionality','color_flag',0,'corrtype','pearson')




%show reciprocity between them


%% Next we were interested in how these things change over time. If things are interacting via low dimensional communication subspace, then should occur at the same time. 
area_name = 'VIS';

[d_time,area_label] = LoadVariable(data,'svca_trialvar_acrosstime',area_name);
rin_time = LoadVariable(data,'rsq_acrosstime',[],15);
ain_time = LoadVariable(data,'auc_acrosstime',[],15);
tavg = LoadVariable(data,'trialavg_timecourse',[]);
tvar = LoadVariable(data,'trialvar_timecourse',[]);
tvar_r = LoadVariable(datar,'trialvar_timecourse_rev',[]);
rout_time = LoadVariable(datar,'rsq_acrosstime_rev',[],15);
aout_time = LoadVariable(datar,'auc_acrosstime_rev',[],15);


%% show example region and motif of activity over time
m = 12;
area_name = 'VIS';
figure; hold on;
Plot_VariableOverTime(squeeze(rin_time(:,m,:)),fp.c_ff,'minmax');
Plot_VariableOverTime(squeeze(tvar(:,m,:)),fp.c_none,'minmax');
ylabel('performance')
yyaxis right
Plot_VariableOverTime(squeeze(ain_time(:,m,:)),fp.c_lr,'minmax');
Plot_VariableOverTime(squeeze(d_time(:,m,:)),fp.c_glm,'minmax');
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])

figure; hold on;
Plot_VariableOverTime(squeeze(rout_time(:,m,:)),fp.c_ff,'minmax');
Plot_VariableOverTime(squeeze(tvar_r(:,m,:)),fp.c_none,'minmax');
ylabel('performance')
yyaxis right
Plot_VariableOverTime(squeeze(aout_time(:,m,:)),fp.c_lr,'minmax');
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])


lm = arrayfun(@(n) fitlm(squeeze(tvar(n,m,:)),squeeze(rin_time(n,m,:))),1:6,'UniformOutput',0);
lm = cellfun(@(x) x.Residuals.Raw,lm,'UniformOutput',0);
rr = cat(2,lm{:})';

lm = arrayfun(@(n) fitlm(squeeze(tvar(n,m,:)),squeeze(ain_time(n,m,:))),1:6,'UniformOutput',0);
lm = cellfun(@(x) x.Residuals.Raw,lm,'UniformOutput',0);
aa = cat(2,lm{:})';

lm = arrayfun(@(n) fitlm(squeeze(tvar(n,m,:)),squeeze(d_time(n,m,:))),1:6,'UniformOutput',0);
lm = cellfun(@(x) x.Residuals.Raw,lm,'UniformOutput',0);
bb = cat(2,lm{:})';

figure; hold on;
Plot_VariableOverTime(rr,fp.c_ff,'minmax');
Plot_VariableOverTime(aa,fp.c_ff,'minmax');
Plot_VariableOverTime(bb,fp.c_glm,'minmax');


t = arrayfun(@(n) corr(rr(n,:)',aa(n,:)','tail','left'), 1:size(aa,1),'UniformOutput',1)

t = arrayfun(@(n) corr(rr(n,:)',aa(n,:)','tail','left'), 1:size(aa,1),'UniformOutput',1)












%% next we were interested what might be driving this subspace
%focus just on the motifs actively engaging that region
% trial_var = LoadVariable(data,'ttt_variability',area_name);
[temp,idx] = sort(nanmean(trial_var),'descend');
x = 1:14;
figure; bar((temp-min(temp))./max(temp-min(temp)));
set(gca,'XTick',x,'XTickLabel',idx)
Z = linkage([nanmean(trial_var)',nanmean(trial_activity)'],'ward');
dendrogram(Z)
trueMotifs = cluster(Z,'MaxClust',2)==1;




%%
%Is this increased with reciprocity?
%scale to the same size
a = squeeze(rin_time(:,m,:));
b = squeeze(rout_time(:,m,:));
routc = b; 
for i = 1:size(a,1)
    coef = polyfit(b(i,:),a(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    routc(i,:) = correction(b(i,:));
end

r = a+routc;

a = squeeze(ain_time(:,m,:));
b = squeeze(aout_time(:,m,:));
aoutc = b; 
for i = 1:size(a,1)
    coef = polyfit(b(i,:),a(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aoutc(i,:) = correction(b(i,:));
end

d = a + aoutc;

a = squeeze(ain_time(:,m,:))./squeeze(tvar(:,m,:));
b = squeeze(aout_time(:,m,:))./squeeze(tvar_r(:,m,:));
aoutc_norm = b; 
for i = 1:size(a,1)
    coef = polyfit(b(i,:),a(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aoutc_norm(i,:) = correction(b(i,:));
end
d_norm = a + aoutc_norm;

%other direction
figure; hold on;
Plot_VariableOverTime(r,fp.c_ff,'minmax');
ylabel('performance')
Plot_VariableOverTime(d,fp.c_lr,'minmax');
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])

%Example plot over time of subspace dimensionality and performance
figure; hold on;
y = Plot_VariableOverTime(r,fp.c_ff,'minmax');
ylabel('performance') 
yyaxis right
x = Plot_VariableOverTime(d_norm,fp.c_glm,'minmax');
ylabel('dimensionality'); ylim([0 1])
set(gca,'YColor','k')
fp.FormatAxes(gca); box on
fp.FigureSizing(gcf,[3 2 6 6],[10 10 14 10])

%get the correlation 
arrayfun(@(n) corr(r(n,:)',d_norm(n,:)','tail','left'), 1:size(x,1),'UniformOutput',1)









%% reciprocity
m = 5;
r = [];
d = [];
d_loc = [];
dnorm = [];
d_loc_norm=[];
for m = 1:14
%get the in and out performance to match scale
a = squeeze(rin_time(:,m,:));
b = squeeze(rout_time(:,m,:));
routc = b; 
for i = 1:size(a,1)
    coef = polyfit(b(i,:),a(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    routc(i,:) = correction(b(i,:));
end

a = squeeze(ain_time(:,m,:));
b = squeeze(aout_time(:,m,:));
aoutc = b; 
for i = 1:size(a,1)
    coef = polyfit(b(i,:),a(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aoutc(i,:) = correction(b(i,:));
end

temp = (routc+squeeze(rin_time(:,m,:))).*squeeze(tvar(:,m,:));
r(:,:,m) = temp./max(temp,[],2);

temp = (squeeze(ain_time(:,m,:))+aoutc)./squeeze(d_time(:,m,:));
dnorm(:,:,m) = temp./max(temp,[],2);

temp = squeeze(ain_time(:,m,:))+aoutc;
d(:,:,m) = temp./max(temp,[],2);

temp = squeeze(d_time(:,m,:));
d_loc(:,:,m) = temp./max(temp,[],2);

end
% 
% r = nanmean(r,3);
% d = nanmean(d,3);
% d_loc = nanmean(d_loc,3);
% dnorm = nanmean(dnorm,3);

m=3;
r = squeeze(r(:,:,m));
d = squeeze(d(:,:,m));
d_loc = squeeze(d_loc(:,:,m));

%basically, when engaged you enter this period of higher dimensionality... 
%but this increases

%%
figure; hold on;
y=r;
x = 1:size(y,2);
errorbar(x,nanmean(y),sem(y))

y=d_loc;
errorbar(x,nanmean(y),sem(y))

yyaxis right
y=d;
x = 1:size(y,2);
errorbar(x,nanmean(y),sem(y))







%%
%plot magnitude of interaction by subspace
figure; hold on;
y=r;
y = y./max(y,[],2);
x = 1:size(y,2);
plot(x,nanmean(y))

%plot the local dimensionality
y = squeeze(d_time(:,m,:));
y = y./max(y,[],2);
x = 1:size(y,2);
plot(x,nanmean(y))


%plot magnitude of interaction by subspace
yyaxis right
y=d;
y = y./max(y,[],2);
x = 1:size(y,2);
plot(x,nanmean(y))



%first do for in prediction | may need to normalize by the trial
%variance... or use a specific dimension. 


%plot the over time plots (try without any weighting or normalization
%plots) ... ensures a consistent subspace over time (and the same order). 

%interestingly, this was also the case with the inverse (

%%THEN SHOW THE EFFECT OF reciprocity on this subspace

%then do the local dimensionality... show that it doesn't really have much
%going on during the subspace period


%% 
Plot_CorrelateValuesBetweenRecordings(aucin,squeeze(localauc(:,:,8)),'average',fp,'right','xlabel','dimensionality','ylabel','performance','color_flag',0)
Plot_CorrelateValuesBetweenRecordings(pin,squeeze(localauc(:,:,8)),'average',fp,'right','xlabel','dimensionality','ylabel','performance','color_flag',0)

%% Next we were interested in how the dimensionality of these may relate to how extensively these regions were interacting. 
area_name = 'VIS';
localdim = LoadVariable(data,'svca',[],0.05);
%correlation between dim and rsq (with an wihtout trial variance)
aucin = LoadVariable(data,'rrr_auc',area_name,15);
trial_var = LoadVariable(data,'ttt_variability',area_name);
pin = LoadVariable(data,'ridge_performance',area_name);
% pin = pin.*trial_var;
%can include per recording
Plot_CorrelateValuesBetweenRecordings(aucin,pin,'average',fp,'left','xlabel','dimensionality','ylabel','performance','color_flag',0)

%same for the recipricol model
aucout = LoadVariable(datar,'rrr_auc',area_name,15);
r_trial_var = LoadVariable(datar,'ttt_variability',area_name);
pout = LoadVariable(datar,'ridge_performance',area_name);
% pout = pout.*trial_var;
Plot_CorrelateValuesBetweenRecordings(aucout,pout,'average',fp,'left','xlabel','dimensionality','ylabel','performance','color_flag',0)

%Test wether it's the reciprocity that matters
poutc = pout; 
for i = 1:size(pin,1)
    coef = polyfit(pout(i,:),pin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    poutc(i,:) = correction(pout(i,:));
end
aucoutc = aucout; 
for i = 1:size(pin,1)
    coef = polyfit(aucout(i,:),aucin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aucoutc(i,:) = correction(aucout(i,:));
end

r = (poutc+pin).*trial_var;
d = (aucoutc+aucin)./nanmean(localdim,3);
Plot_CorrelateValuesBetweenRecordings(d,r,'average',fp,'left','xlabel','dimensionality','ylabel','performance','color_flag',0)











%% Get the local dimensionality over time 
area_name = 'RSP';
m=5;
[locald,area_label] = LoadVariable(data,'svca_trialvar_acrosstime',[]);
[rsq,area_label] = LoadVariable(data,'rsq_acrosstime',[],30);
[rsqout,area_label] = LoadVariable(datar,'rsq_acrosstime_rev',[],30);
[aucin,area_label] = LoadVariable(data,'auc_acrosstime',[],30);
[aucout,area_label] = LoadVariable(datar,'auc_acrosstime_rev',[],30);
tavg = LoadVariable(data,'trialavg_timecourse',[]);

%%
area_name = 'VIS';
% close all;
figure; hold on; 
t = tiledlayout(4,4); t.TileSpacing = 'compact'; t.Padding = 'compact';
for i = 1:14
    nexttile; hold on;
    x = squeeze(locald(:,i,strcmp(area_label,area_name)==1,:));
    x = x./nanmax(x,[],2);
    plot(nanmean(x))    
    y = squeeze(rsq(:,i,strcmp(area_label,area_name)==1,:));
    y = y./nanmax(y,[],2);
    plot(nanmean(y),'r')    
    y = squeeze(rsqout(:,i,strcmp(area_label,area_name)==1,:));
    y = y./nanmax(y,[],2);
    plot(nanmean(y),'b')     
    plot(nanmean(y),'r')    
    y = squeeze(tavg(:,i,strcmp(area_label,area_name)==1,:));
    y = y./nanmax(y,[],2);
    plot(nanmean(y),'k')  
   yyaxis right
    y = squeeze(aucin(:,i,strcmp(area_label,area_name)==1,:));
    y = y./nanmax(y,[],2);
    plot(nanmean(y),'c')    
    y = squeeze(aucout(:,i,strcmp(area_label,area_name)==1,:));
    y = y./nanmax(y,[],2);
    plot(nanmean(y),'g')   
end

%% correlations to test
% per recording and per motif 
abslag=[];
truelag=[];
for i = 1:8
    verbose=0;
    area_name = area_label{i};
    a = squeeze(locald(:,:,strcmp(area_label,area_name)==1,:));
    b = squeeze(rsq(:,:,strcmp(area_label,area_name)==1,:));
    c = squeeze(rsqout(:,:,strcmp(area_label,area_name)==1,:));
    d = squeeze(tavg(:,:,strcmp(area_label,area_name)==1,:));
    e = squeeze(aucin(:,:,strcmp(area_label,area_name)==1,:));
    f = squeeze(aucout(:,:,strcmp(area_label,area_name)==1,:));

    %cross correlation between in and out performance, per motif
    [rho,lag] = PerformanceXCorr(b,c,1,verbose);
    nanmean(rho(:))
    nanmedian(abs(lag(:)))

    %cross correlation between in performance and aucin
%     [rho,lag] = PerformanceXCorr(e,f,1,1);
%     nanmean(rho(:))
%     nanmedian(abs(lag(:)))

    %cross correlation between trial average and in performance
    [rho,lag] = PerformanceXCorr(d,b,1,verbose);
    nanmean(rho(:))
    nanmedian(abs(lag(:)))


    %cross correlation between trial average and out performance
    [rho,lag] = PerformanceXCorr(d,c,1,verbose);
    nanmean(rho(:))
    nanmedian(abs(lag(:)))


    %cross correlation between avg performance and local dimensionality
    %negative means that it comes after, positive means before. 
    [rho,lag,brho] = PerformanceXCorr((c+b)./2,a,1,verbose);
    nanmean(rho(:))
%     abslag(i,:) = nanmean(brho);
%     truelag(i,:) = nanmean(abs(lag));

    %cross correlation between aucin and local dimensionality
    [rho,lag,brho] = PerformanceXCorr(e,a,1,0);
    nanmean(rho(:))
    abslag(i,:) = nanmean(brho);
    truelag(i,:) = nanmean((lag));

    %the question becomes whether interareal dimensionality and and local
    %dimensionality are correlated... if they are correlated then they both
    %drop when the place is interacting and then open up afterwards. 
       
    
end

%%

ttest(abslag',0)
[p,~,stats] = anova1(truelag')
t = multcompare(stats)

%% for a given brain region, correlate the magnitude of lag with the reciprocity
pin = LoadVariable(data,'ridge_performance',area_name);
pout = LoadVariable(datar,'ridge_performance',area_name);
poutc = pout; 
for i = 1:size(pin,1)
    coef = polyfit(pout(i,:),pin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    poutc(i,:) = correction(pout(i,:));
end

r = (pin+poutc).*nanmean(ntrial,3);


%% The problem is that we have motifs that are always lower dimensional and always less reciprocal across all areas... 
%and this makes no sense because that motif may not even engage that
%region... suggests that there is some other variable (ntrial?) magnitude
%of trial to trial variability? magnitude of engagement?
%is trial-to-trial variability lower for some motifs across all areas? 
ntrial = LoadVariable(data,'ttt_variability',[]);
%zscore across motifs
temp = zscore(ntrial,0,2);

%if topography, then some would not be more than zero or not be correlated








%%
close all
fp = fig_params_deconvolutionpaper;
area_name = 'VIS';
pin = LoadVariable(data,'ridge_performance',area_name);
pout = LoadVariable(datar,'ridge_performance',area_name);
n = LoadVariable(data,'ntrial',area_name);
aucin = LoadVariable(data,'rrr_auc',area_name);
[aucout,area_label] = LoadVariable(datar,'rrr_auc',area_name);
ttt_var = squeeze(ntrial(:,:,strcmp(area_label,area_name)));
ttt_all = nanmean(ntrial,3);
ttt_out = nanmean(ntrial(:,:,strcmp(area_label,area_name)==0),3);
% localdim = LoadVariable(data,'svca_trialvar',[],0.05);
locald = squeeze(localdim(:,:,strcmp(area_label,area_name)));
otherd = nanmean(localdim(:,:,strcmp(area_label,area_name)==0),3);


%for each recording, get the regression between values so that you can
%scale pin to match pout
poutc = pout; 
for i = 1:size(pin,1)
    coef = polyfit(pout(i,:),pin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    poutc(i,:) = correction(pout(i,:));
end

aucoutc = aucout; 
for i = 1:size(pin,1)
    coef = polyfit(aucout(i,:),aucin(i,:),1); 
    correction = @(x) coef(1)*x+coef(2);
    aucoutc(i,:) = correction(aucout(i,:));
end


%correlated the order of motif dimensnionality across areas
%would be interesting if some motifs were always high dim (like the
%traveling wave) vs more static ones. 

m = [5,8,6,11];
% m = [12, 5, 8]; %ss, vis, and wave
% r = zscore(pin./ttt_var,0,2)+zscore(pout./ttt_var,0,2);
% d = zscore(aucin./otherd,0,2)+zscore(aucout./locald,0,2);

r = abs(pin-poutc); %this takes care of any normalization across recordings because it put everythin  in the same scale
r = r/max(r(:));
d = (aucin+aucoutc)  ./nanmean(localdim,3); %aucin./otherd+aucout./locald;
% d = nanmean(localdim,3); %aucin./otherd+aucout./locald;
r = (pin+poutc).*ttt_all;
x = zscore(r,0,2);
y = zscore(d,0,2); 
x = (r);
% y = (d); 
figure; hold on;
%I think we need to normalize across motifs in a recording
plot(x(:),y(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(x(:,m(n)),y(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('interactivity');
ylabel('relative dimensionality');
[rho,p] = corr(x(:),y(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');

%% figure (plot the marker)
figure; hold on;
x = nanmean(pin.*ttt_all); 
y = nanmean(aucin./nanmean(localdim,3)); 
%I think we need to normalize across motifs in a recording
plot(x(:),y(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(x(:,m(n)),y(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('interactivity');
ylabel('relative dimensionality');
[rho,p] = corr(x(:),y(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');


figure; hold on;
x = nanmean(poutc.*ttt_all); 
y = nanmean(aucout./nanmean(localdim,3)); 
%I think we need to normalize across motifs in a recording
plot(x(:),y(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(x(:,m(n)),y(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('interactivity');
ylabel('relative dimensionality');
[rho,p] = corr(x(:),y(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');



%% Confirm that motifs are ordered differently dimensionality and performance across areas (i.e. not some global effect coming in). 
a = []; b= [];
for j = 1:8
    aucin = LoadVariable(data,'rrr_auc',area_label{j});
    aucout = LoadVariable(datar,'rrr_auc',area_label{j});    
    pin = LoadVariable(data,'ridge_performance',area_label{j});
    pout = LoadVariable(datar,'ridge_performance',area_label{j});
    ttt_var = squeeze(ntrial(:,:,strcmp(area_label,area_label{j})));
    ttt_all = nanmean(ntrial,3);
    ttt_out = nanmean(ntrial(:,:,strcmp(area_label,area_label{j})==0),3);


    %for each recording, get the regression between values so that you can
    %scale pin to match pout
    poutc = pout; 
    for i = 1:size(pin,1)
        coef = polyfit(pout(i,:),pin(i,:),1); 
        correction = @(x) coef(1)*x+coef(2);
        poutc(i,:) = correction(pout(i,:));
    end

    aucoutc = aucout; 
    for i = 1:size(pin,1)
        coef = polyfit(aucout(i,:),aucin(i,:),1); 
        correction = @(x) coef(1)*x+coef(2);
        aucoutc(i,:) = correction(aucout(i,:));
    end
    d = (aucin+aucoutc);%./nanmean(localdim,3); %aucin./otherd+aucout./locald;
    r = (pin+poutc);%.*ttt_all;    
    
    [~,a(j,:)] = sort(nanmean(zscore(r,0,2)),'ascend'); %increasing interactivity
    [~,b(j,:)] = sort(nanmean(zscore(d,0,2)),'ascend'); %increasing dimensionality 
    
end
temp = tril(corr(a',a','type','spearman'),-1);
temp(temp==0)=NaN;
figure; hold on; 
histogram(temp(:),'BinWidth',0.1);

temp = tril(corr(b',b','type','spearman'),-1);
temp(temp==0)=NaN;
figure; hold on; 
histogram(temp(:),'BinWidth',0.1);

% also worth highlighting the pretty considerable variability in dimensionality you have on a local level across motifs. 


%% what are cross population things that are playing a role?
% one potential is reciprocol interactions? 
%define measure of reciprocity = pev in * pev out
close all
fp = fig_params_deconvolutionpaper;
area_name = 'VIS';
pin = LoadVariable(data,'ridge_performance',area_name);
pout = LoadVariable(datar,'ridge_performance',area_name);
ntrial = LoadVariable(data,'ntrial',area_name);
aucin = LoadVariable(data,'rrr_auc',area_name);
aucout = LoadVariable(datar,'rrr_auc',area_name);
localdim = LoadVariable(data,'svca_trialvar',area_name,0.05);
m = [5,8,11];
% m = [12, 5, 8]; %ss, vis, and wave
d = zscore(aucin,0,2)+zscore(aucout,0,2);
r = ttt_var./pin+ttt_var./pout;


%%
figure; hold on;
%I think we need to normalize across motifs in a recording
plot(r(:),d(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(r(:,m(n)),d(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('reciprocity');
ylabel('dimensionality');
[rho,p] = corr(r(:),d(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');

% pinc = pin.*ttt_var;
% poutc = poutc.*ttt_out; 
% for i = 1:size(pin,1)
%     coef = polyfit(poutc(i,:),pin(i,:),1); 
%     correction = @(x) coef(1)*x+coef(2);
%     poutc(i,:) = correction(poutc(i,:));
% end

%I don't think we can compare the rsq between motifs. Within a motif, we
%can compare the reciprocity. 
% r = (pinc+poutc).*ttt_all;
r = (pinc+poutc).*ttt_all;

%compare predictability in and out for three motifs...
figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(ttt_var,m,fp);
ylabel('variance in');

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(ttt_out,m,fp);
ylabel('variance out');

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(pin,m,fp);
ylabel('r^2 In');

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(poutc,m,fp);
ylabel('r^2 Out');

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(pin.*ttt_var,m,fp);
ylabel('r^2 In');

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(poutc.*ttt_out,m,fp);
ylabel('r^2 Out');
% set(gca,'ylim',[0.01 0.06])

%and their resulting reciprocity...
figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(r,m,fp);
ylabel('Interactivity');
% set(gca,'ylim',[0.01 0.06])

%and their in and out dimensionality
figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(aucin,m,fp);
ylabel('AUC In');
set(gca,'ylim',[0.8 1])

figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(aucout,m,fp);
ylabel('AUC Out');
set(gca,'ylim',[0.8 1])

%compare local dimensionality
figure; hold on; 
fh = Plot_CompareValueBetweenMotifs(localdim,m,fp);
ylabel('Local Dim');

% find the proper normalization method (if any)
x = nanmean(pin);
y = nanmean(pout);
x = pin./max(pin,[],2);
y = pout./max(pout,[],2);
x = nanmean(zscore(pin,0,2));
y = nanmean(zscore(pout,0,2));
x = zscore(pin,0,2)+zscore(pout,0,2); 
y = zscore(pout,0,2);
x = zscore(aucin,0,2)+zscore(aucout,0,2);
y =  zscore(aucout,0,2);
figure; hold on;
%I think we need to normalize across motifs in a recording
plot(x(:),y(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(x(:,m(n)),y(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('PEV in');
ylabel('PEV out');
[rho,p] = corr(x(:),y(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');




figure; hold on;
%I think we need to normalize across motifs in a recording
plot(r(:),d(:),'marker','o','linestyle','none');
%plot our above motifs on there
col = getColorPalet(numel(m));
p = arrayfun(@(n) plot(r(:,m(n)),d(:,m(n)),'marker','o','linestyle','none','markerfacecolor',col(n,:)), 1:numel(m),'UniformOutput',0);
temp = arrayfun(@(x) num2str(x),m,'UniformOutput',0);
legend([p{:}],temp);
xlabel('reciprocity');
ylabel('dimensionality');
[rho,p] = corr(r(:),d(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');


%mean version of the chart should be similar
figure; hold on;
r = abs(pin-pout).*(pin+pout);
%I think we need to normalize across motifs in a recording
plot(nanmean(r),nanmean(d),'marker','o','linestyle','none');
%plot our above motifs on there
xlabel('reciprocity');
ylabel('dimensionality');
[rho,p] = corr(nanmean(r)',nanmean(d)','type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');


%compare subspace dimensionality and local dimensionality
figure; hold on;
plot(d(:),localdim(:),'marker','o','linestyle','none');
xlabel('subspace dimensionality');
ylabel('local dimensionality');
[rho,p] = corr(d(:),localdim(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');

%compare reciprocity local dimensionality
figure; hold on;
plot(r(:),localdim(:),'marker','o','linestyle','none');
xlabel('reciprocity');
ylabel('local dimensionality');
[rho,p] = corr(r(:),localdim(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_name,rho,p),'fontweight','normal');

corr(aucin(:),ntrial(:))
corr(aucout(:),ntrial(:))
[rho,p]=partialcorr(r(:),d(:),ntrial(:))

%% reciprocal predictions also low dimensional? 
figure; hold on; 
ndim = 15;
thresh = 0.80;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
rrr_mdl = LoadVariable(data,'rel_performance','SS');
rrr_mdl = squeeze(rrr_mdl(:,m,1:15));
arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(rrr_mdl,1));
plot(1:size(rrr_mdl,2), nanmean(rrr_mdl),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);  
d = arrayfun(@(n) find(rrr_mdl(n,:)>=thresh,1,'first'),1:size(rrr_mdl,1));
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel('performance');
title(sprintf('Reverse | Motif %d needed %d to %d dim to capture %0.1d%% PEV',m,min(d), max(d),thresh*100),'Fontweight','normal')

%% Expand to other brain areas... it changes!! i.e. topography... why? 
%one hypothesis is that this is due to the structure of the local circuit?

%look at beta weights in 'prelimbic' area































%% does dimensionality scale with engagement
%get AUC of the curves
area_name = 'PRE';
auc = LoadVariable(data,'rrr_auc',area_name);
aucr = LoadVariable(datar,'rrr_auc',area_name);
engment = LoadVariable(data,'engagement',area_name);
engment_norm = zscore(engment,0,2); %(engment)./(max(engment,[],2)); %normalize by the max prediction and activity across motifs per recording

close all; figure('units','normalized','position',[0.1 0.3 0.6 0.3]); hold on;
t=tiledlayout(1,3); t.TileSpacing = 'compact'; t.Padding = 'compact';
nexttile
rrr_d = LoadVariable(data,'rrr_dim',area_name);
rrr_d_reverse = LoadVariable(datar,'rrr_dim',area_name);    
x = rrr_d(:);
y = rrr_d_reverse(:);
x = x+(rand(numel(x),1)/4);
y = y+(rand(numel(y),1)/4);
plot(x(:),y(:),'o','markersize',3,'color',[0.8 0.1 0.1])
[rho,p] = corr(rrr_d(:),rrr_d_reverse(:),'type','spearman','rows','complete');
xlabel('# Dim to Reach 80% PEV Incoming');
ylabel('# Dim to Reach 80% PEV Outgoing');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');

nexttile
plot(auc(:),aucr(:),'marker','o','linestyle','none');
xlabel('AUC Incoming');
ylabel('AUC Outgoing');
[rho,p] = corr(auc(:),aucr(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');

pev = LoadVariable(data,'ridge_performance',area_name);
pevr = LoadVariable(datar,'ridge_performance',area_name);
nexttile
plot(pev(:),pevr(:),'marker','o','linestyle','none');
xlabel('PEV (full model) Incoming');
ylabel('PEV (full model) Outgoing');
[rho,p] = corr(pev(:),pevr(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');




auc_norm = zscore(auc,0,2);
aucr_norm = zscore(aucr,0,2);

close all; figure('units','normalized','position',[0.1 0.3 0.6 0.3]); hold on;
t=tiledlayout(1,3); t.TileSpacing = 'compact'; t.Padding = 'compact';
nexttile
rrr_d = LoadVariable(data,'rrr_dim',area_name);
rrr_d_reverse = LoadVariable(datar,'rrr_dim',area_name);    
x = rrr_d(:);
y = rrr_d_reverse(:);
x = x+(rand(numel(x),1)/4);
y = y+(rand(numel(y),1)/4);
plot(x(:),y(:),'o','markersize',3,'color',[0.8 0.1 0.1])
[rho,p] = corr(rrr_d(:),rrr_d_reverse(:),'type','spearman','rows','complete');
xlabel('# Dim to Reach 80% PEV Incoming');
ylabel('# Dim to Reach 80% PEV Outgoing');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');

nexttile
plot(auc_norm(:),aucr_norm(:),'marker','o','linestyle','none');
xlabel('AUC Incoming');
ylabel('AUC Outgoing');
[rho,p] = corr(auc_norm(:),aucr_norm(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');

pev = LoadVariable(data,'ridge_performance',area_name);
pevr = LoadVariable(datar,'ridge_performance',area_name);
nexttile
plot(pev(:),pevr(:),'marker','o','linestyle','none');
xlabel('PEV (full model) Incoming');
ylabel('PEV (full model) Outgoing');
[rho,p] = corr(pev(:),pevr(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');


x = [auc(:),aucr(:)];
y = engment(:);
x = [auc_norm(:),aucr_norm(:)];
y = engment_norm(:);
lm = fitlm(x,y,'interactions')
% CVMdl = crossval(lm) 

x = [auc(:),aucr(:),y(:)]';
xhat = mean(x,2);
A = x-xhat;
[U,S,~] = svd(A);
d = U(:,1);
t = d'*A;
t1 = min(t);
t2 = max(t);
xhat = xhat + [t1,t2].*d; % size 3x2
xhat = xhat';

figure; hold on; 
plot3(auc(:),aucr(:),engment(:),'Marker','o','linestyle','none','color',[0.8 0.1 0.1]);
% plot3(xhat(:,1),xhat(:,2),xhat(:,3),'Marker','none','linestyle','-','color','k','linewidth',2);
xlabel('auc in'); ylabel('auc out'); zlabel('engagement');

figure; hold on;
reciprocity = auc_norm(:)+aucr_norm(:);
plot(reciprocity,engment_norm(:),'marker','o','linestyle','none')
[rho,p] = corr(reciprocity,engment_norm(:),'type','pearson','rows','complete');
title(sprintf('Area %s Rho %0.2f pval %0.2f',area_label{i},rho,p),'fontweight','normal');
xlabel('reciprocity (auc in + auc out)'); ylabel('engagement');


%% Do the beta weights match what one might expect?
%full beta weights for Vis motif for example recording
area_name = 'VIS';
m = [5,8,11];
b = LoadVariable(data,'rrr_beta',area_name,1);
[reg,area_lbl] = LoadVariable(data,'beta_region',area_name);
%full beta weights
figure('units','normalized','position',[0 0 1 1]); hold on; 
t=tiledlayout(4,3); t.TileSpacing = 'compact'; t.Padding = 'compact';
for cur_m = 1:numel(m)
    nexttile; hold on;
    x = squeeze(b(1,m(cur_m),:));
    bar(x,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(1,m(cur_m),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title('raw beta weights','fontweight','normal');
    xlabel('neurons'); ylabel('Beta Weights');
end

%mean fr normalized beta weights
fr = LoadVariable(data,'mean_fr',area_name);
b = LoadVariable(data,'rrr_beta',area_name,1);
for cur_m = 1:numel(m)
    nexttile; hold on;
    x = b.*fr;
    x = squeeze(x(1,m(cur_m),:));
    bar(x,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(1,m(cur_m),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title('Beta x Mean FR','fontweight','normal');
    xlabel('neurons'); ylabel('Weights');
end

%mean fr normalized beta weights
fr = LoadVariable(data,'trialvar_fr',area_name);
b = LoadVariable(data,'rrr_beta',area_name,1);
for cur_m = 1:numel(m)
    nexttile; hold on;
    x = b.*fr;
    x = squeeze(x(1,m(cur_m),:));
    bar(x,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(1,m(cur_m),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title('Beta x |Trial2trial FR|','fontweight','normal');
    xlabel('neurons'); ylabel('Weights');    
end

%mean fr normalized beta weights
b = LoadVariable(data,'rrr_beta_weightedsorted',area_name,1);
for cur_m = 1:numel(m)
    nexttile; hold on;
    x = squeeze(b(1,m(cur_m),:));
    bar(x,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(1,m(cur_m),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title('Beta | weighted and sorted','fontweight','normal');
    xlabel('neurons'); ylabel('Weights');
end

title(t,sprintf('area %s motifs %d, %d, and %d',area_name, m),'fontweight','normal')

%%

%show average across all six recordings
m=11;
[b,area_lbl] = LoadVariable(data,'rrr_bsum_max',area_name,1);
area_lbl(ismember(area_lbl,area_name))=[];
figure; hold on; 
plot([1,numel(area_lbl)],[0,0],'linewidth',1.5,'linestyle',':','color','k')
x = squeeze(b(:,m,:));
arrayfun(@(n) plot(1:size(x,2),x(n,:),'marker','o','markersize',3,'linewidth',0.5,'color',[0.75 0.75 0.75]), 1:size(x,1));
errorbar(1:size(x,2), nanmean(x),sem(x),'linewidth',1.5,'color',[0.80,0.1,0.1],'marker','o','markersize',5);  
set(gca,'xtick',1:numel(area_lbl),'xticklabel',area_lbl,'XTickLabelRotation',45)
ylabel('Summed Normalized Beta Weigths'); xlabel('Source Area');
p = anova1(x,[],'off');
title(sprintf('Contributions to %s prediction p=%0.2f | motif %d',area_name,p,m),'fontweight','normal')

% Higher beta weights scale with more engagement of that region
% the more I engage RSP the more beta (when controlling for the engagement of VIS)
source_region = 'RSP';
[b,area_label] = LoadVariable(data,'rrr_bsum_max',area_name,1);
engment = LoadVariable(data,'engagement',[]);

%across all motifs and recordings
x = squeeze(b(:,:,strcmp(area_label,source_region)==1));
y = squeeze(engment(:,:,strcmp(area_label,source_region)==1));
yy = squeeze(engment(:,:,strcmp(area_label,area_name)==1));
[rho,p]= partialcorr(y(:),x(:),yy(:),'rows','complete','tail','right');
yr = fitlm(yy(:),y(:));
xr = fitlm(yy(:),x(:));
figure; hold on; 
xr = xr.Residuals.Standardized(~isnan(xr.Residuals.Standardized));
yr = yr.Residuals.Standardized(~isnan(yr.Residuals.Standardized));
plot(xr,yr,'marker','o','markersize',3,'linestyle','none','color',[0.8, 0.1 0.1])
P = polyfit(xr,yr,1);
yfit = P(1)*xr+P(2);
plot(xr,yfit,'r-','linewidth',2);
ylabel(sprintf('Predictive Weight of %s area',source_region)); xlabel(sprintf('Engagement of %s area',source_region));
title(sprintf('More engagement = stronger Beta | Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')


%% Cool, so that adds up to a subspace... now let's look
























%% If we focus in on the motif betas... what do they look like?
%plot example betas on two areas and two motifs show choice of normalizations
close all
m = [5,8,10,2];
area_temp = {'VIS','MOs'};
b = LoadVariable(data,'rrr_beta',area_temp{1},1);
a = LoadVariable(data,'rrr_beta',area_temp{2},1);
b(2,:,:)=a(1,:,:); b(3:end,:,:)=[];
[reg,area_lbl] = LoadVariable(data,'beta_region',area_temp{1});
a = LoadVariable(data,'beta_region',area_temp{2});
reg(2,:,:)=a(1,:,:); reg(3:end,:,:)=[];

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,numel(m));
idx = combvec(m,[1,2]); %motifs and recordings
for j = 1:size(idx,2)
    nexttile; hold on;
    bar(squeeze(b(idx(2,j),idx(1,j),:)),'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(idx(2,j),idx(1,j),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title(sprintf('area %s motif %d',area_temp{idx(2,j)},idx(1,j)),'fontweight','normal')
    xlabel('neurons'); ylabel('Beta Weights');
end
title(t,'Example beta weights (dimension 1)')

%normalize betas by strongest
b = LoadVariable(data,'rrr_beta_max',area_temp{1},1);
a = LoadVariable(data,'rrr_beta_max',area_temp{2},1);
b(2,:,:)=a(1,:,:); b(3:end,:,:)=[];
[reg,area_lbl] = LoadVariable(data,'beta_region',area_temp{1});
a = LoadVariable(data,'beta_region',area_temp{2});
reg(2,:,:)=a(1,:,:); reg(3:end,:,:)=[];

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,numel(m));
idx = combvec(m,[1,2]); %motifs and recordings
for j = 1:size(idx,2)
    nexttile; hold on;
    bar(squeeze(b(idx(2,j),idx(1,j),:)),'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5,'EdgeAlpha',0);
    lbl = squeeze(reg(idx(2,j),idx(1,j),:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(find(lbl==lblu(i)),zeros(1,sum(lbl==lblu(i))),'color',c(i,:),'linewidth',2);
    end
    legend([p{:}],area_lbl(lblu));
    title(sprintf('area %s motif %d',area_temp{idx(2,j)},idx(1,j)),'fontweight','normal')
    xlabel('neurons'); ylabel('Normalized Beta Weights');
end
title(t,'Normalizing beta weights to max')


%% Plot the full beta weights | 3 motifs for 3 target areas across 10 dim
% so it appears that some regions, are using the same relationships just
% scaled up/down
m = [5,8,6];
area_name = {'VIS','RSP','MOs'};
b_all = cell(1,10);
for cur_d = 1:10
   [b,area_lbl] = LoadVariable(data,'rrr_beta',[],cur_d);
   b_all{cur_d} = squeeze(nanmean(b(:,m,ismember(area_lbl,area_name),:),1));
end
b_all = cat(4,b_all{:});
close all; 
figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on;
t = tiledlayout(numel(m),numel(area_name)); t.TileSpacing = 'compact'; t.Padding = 'compact';
for cur_m = 1:numel(m)    
    for cur_a = 1:numel(area_name)
        x = squeeze(b_all(cur_m,cur_a,:,:));
        x(sum(isnan(x),2)>0,:)=[];
        nexttile; hold on; 
        imagesc(x,[-0.1 0.1]); colormap redgreencmap; c=colorbar;
        ylabel(c,'Normalized beta');
        title(sprintf('Area: %s | Motif %d',area_name{cur_a},m(cur_m)),'fontweight','normal');
        xlabel('predictive dimensions'); ylabel('Neurons');
        axis square
        xlim([0.5, cur_d+0.5])
        ylim([0.5 size(x,1)+0.5]);
    end        
end

%% plot the relative engagement of each 
%do average beta weights match what we might expect? 
% close all;
area_name = 'VIS';
m = [5,7,8,6,2];
for cur_d = 1:5
    [b,area_lbl] = LoadVariable(data,'rrr_bsum_max',area_name,cur_d);
    area_lbl(strcmp(area_lbl,area_name))=[];
    %choose three motifs that all engage visual area    
    figure('position',[ 229 559 1538 420]); hold on; 
    t = tiledlayout(1,numel(m)); t.TileSpacing = 'compact'; t.Padding = 'compact';
    yval = [floor(nanmin(b(:,m,:),[],'all')),ceil(nanmax(b(:,m,:),[],'all'))];
    for i = 1:numel(m)
       nexttile; hold on
       plot([1,numel(area_lbl)],[0,0],'linewidth',1.5,'linestyle',':','color','k')
       x = squeeze(b(:,m(i),:));
       arrayfun(@(n) plot(1:size(x,2),x(n,:),'marker','o','markersize',3,'linewidth',0.5,'color',[0.75 0.75 0.75]), 1:size(x,1));
       errorbar(1:size(x,2), nanmean(x),sem(x),'linewidth',1.5,'color',[0.80,0.1,0.1],'marker','o','markersize',5);  
       set(gca,'ylim',yval,'xtick',[1:numel(area_lbl)],'xticklabel',area_lbl,'XTickLabelRotation',45)
       ylabel('Summed Normalized Beta Weigths');
       xlabel('Source Area');
       title(sprintf('Motif %d',m(i)),'FontWeight','normal');
    end
    title(t,sprintf('Source area contributions to %s | dimension %d',area_name,cur_d),'fontweight','normal')
end %dimension loop 

%% Make our life easy and plot 3 motifs and 3 areas across 10 dim
m = [5,8,6];
area_name = 'VIS';
source_areas = {'MOs','RSP','SS'};
b_all = cell(1,10);
for cur_d = 1:10
   [b,area_lbl] = LoadVariable(data,'rrr_bsum_max',area_name,cur_d);
   area_lbl(strcmp(area_lbl,area_name))=[];
   source_idx = ismember(area_lbl,source_areas);
   b_all{cur_d} = b(:,m,source_idx);
end
b_all = cat(4,b_all{:});

%per source, plot the three motifs across dimensions
close all
figure('position',[208 559 1403 420]); hold on; 
t = tiledlayout(1,numel(source_areas)); t.TileSpacing = 'compact'; t.Padding = 'compact';
c = getColorPalet(numel(m));
for cur_a = 1:numel(source_areas)
    nexttile; hold on;
    p = cell(1,numel(m));
    lbl = cell(1,numel(m));
    for cur_m = 1:numel(m)
        x = squeeze(b_all(:,cur_m,cur_a,:));
        shadedErrorBar(1:size(x,2),nanmean(x),sem(x),'lineprops',{'color',c(cur_m,:)});
        p{cur_m} = plot(1:size(x,2), nanmean(x), 'linewidth',1.5,'color',c(cur_m,:));
        lbl{cur_m} = sprintf('Motif %d',m(cur_m));
    end
    pval = arrayfun(@(n) anova1(squeeze(b_all(:,:,cur_a,n)),[],'off'), 1:size(b_all,4))    
    title(sprintf('Source area: %s',source_areas{cur_a}),'fontweight','normal');
    xlabel('predictive dimensions'); ylabel('Summed normalized beta');
    legend([p{:}],lbl);
end


%% How similar are betas across dimensions per area
%error bar of similarity between subsequent dimensions
close all
m = [5,8,6];
temp = NaN(8,14,10,6);
for cur_d = 1:10
   [b,area_lbl] = LoadVariable(data,'rrr_beta_max',[],cur_d);
   [b2,area_lbl] = LoadVariable(data,'rrr_beta_max',[],cur_d+1);
   %for each motif and each area, get the similarity in subsequent dimesnions
   for cur_a = 1:size(b,3)
       for cur_m = 1:size(b,2)
            x = squeeze(b(:,cur_m,cur_a,:));
            bad_idx = sum(~isnan(x),2)==0;
            x(bad_idx,:)=[];
            y = squeeze(b2(:,cur_m,cur_a,:));
            y(bad_idx,:)=[];
            rho = corr(x',y','rows','complete');
            temp(cur_a,cur_m,cur_d,bad_idx==0) = diag(rho);
       end
   end
end
figure('units','normalized','position',[0.2 0.2 0.6 0.6]); hold on; 
t = tiledlayout(4,2); t.TileSpacing = 'compact'; t.Padding = 'compact';
for cur_a = 1:size(temp,1)
    nexttile; hold on; 
    for cur_m = 1:size(temp,2) 
        x = squeeze(temp(cur_a,cur_m,:,:));
%         errorbar(1:size(x,1),nanmean(x,2),sem(x,2),'marker','o','markersize',5,'color',[0.8 0.1 0.1]);
        plot(1:size(x,1),nanmean(x,2),'marker','none','linewidth',1,'color',[0.8 0.1 0.1]);         
    end
%     legend()
    title(sprintf('Area %s',area_lbl{cur_a}));
end


%

%% How similar are betas across our motifs per area



















end %function end










