function AnalyzeCCs_draft2()
% camden - timeless

%Load the data_rrr for our specific motifs from each rec
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';

data = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPED\w*.mat'],0,{folder}); 
    %get motif '6,7,8'
%     idx = contains(fn,'motif6') | contains(fn,'motif7') | contains(fn,'motif8');
    data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif'),fn);
end

datacca = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*CCA_muaflag1_GROUPED\w*.mat'],0,{folder}); 
    datacca{cur_rec} = cellfun(@(x) load(x,'r','area_label','pairings','motif'),fn);
end

datar = cell(1,6); %reversed
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDREVERSE\w*.mat'],0,{folder}); 
    datar{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','area_val','paired_areas','rrr_V','rrr_B','grouping','motif'),fn);
end


%remove noise motif and null motif
for i = 1:6
%    idx = ismember(arrayfun(@(n) datacca{i}(n).motif, 1:size(datacca{i},2)),[2,16]);
%    datacca{i}(idx)=[];   
   idx = ismember(arrayfun(@(n) datar{i}(n).motif, 1:size(datar{i},2)),[2,16]);
   datar{i}(idx)=[];      
%    idx = ismember(arrayfun(@(n) data{i}(n).motif, 1:size(data{i},2)),[2,16]);
%    data{i}(idx)=[];     
end


%% Question 1: When brain area is engaged by a motif, does it have more predictable varaibilty with other areas in the brian? 
m = 5;
mout =6;
area_name = 'VIS';
full_mdl = LoadVariable(data,'ridge_performance',area_name);

%plot example recording Mean +/- SEM for full model 
figure; hold on; 
errorbar(1, nanmean(full_mdl(:,m)),sem(full_mdl(:,m)),'linestyle','none','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
plot(cat(1,ones(1,size(full_mdl,1)),2*ones(1,size(full_mdl,1))),full_mdl(:,[m,mout])','marker','o','markersize',5,'color',[0.65 0.65 0.65],'linestyle','-')
errorbar(2,nanmean(full_mdl(:,mout)),sem(full_mdl(:,mout)),'linestyle','none','marker','o','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
set(gca,'xlim',[0.5 2.5],'ylim',[0.07 0.13],'xtick',[1,2],'xticklabel',[m,mout])
xlabel('motif'); ylabel('performance (10fold xval)');
title(sprintf('Example Performance motif %d area %s',m,area_name),'fontweight','normal')

%'vis' motifs should have more activity in it than, say 'non-visual' motif X
engment = LoadVariable(data,'engagement',area_name);
figure; hold on; 
bar(1,nanmean(engment(:,m),'all'),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0)
bar(2,nanmean(engment(:,mout),'all'),'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8],'EdgeAlpha',0)
for i = 1:size(engment,1)
   plot([1,2],[nanmean(engment(i,m)),nanmean(engment(i,mout))],'color',[0.25 0.25 0.25],'marker','o','markersize',5)
end
[~,p] = ttest(nanmean(engment(:,m),2),nanmean(engment(:,mout),2),'tail','right');
ylabel('Activity (normalized FR)')
xlabel(sprintf('motif p=%0.2f',p));
title(sprintf('Average "Engagement" of %s area',area_name),'fontweight','normal')
set(gca,'xtick',[1,2],'xticklabel',[m,mout]);

%'vis' motifs should have better performance than, say 'non-visual' motif X
figure; hold on; 
bar(1,nanmean(full_mdl(:,m),'all'),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0)
bar(2,nanmean(full_mdl(:,mout),'all'),'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8],'EdgeAlpha',0)
for i = 1:size(full_mdl,1)
   plot([1,2],[nanmean(full_mdl(i,m)),nanmean(full_mdl(i,mout))],'color',[0.25 0.25 0.25],'marker','o','markersize',5)
end
[~,p] = ttest(nanmean(full_mdl(:,m),2),nanmean(full_mdl(:,mout),2),'tail','right');
ylabel('performance')
xlabel(sprintf('motif p=%0.2f',p));
title('Predicting VIS trial-trial activity from other areas','fontweight','normal')
set(gca,'xtick',1:2,'xticklabel',[m,mout]);

%to combine across all recordings, normalize engagement and prediction
engment_norm = (engment-min(engment,[],2))./(max(engment,[],2)-min(engment,[],2)); %normalize by the max prediction and activity across motifs per recording
full_mdl_norm =  (full_mdl-min(full_mdl,[],2))./(max(full_mdl,[],2)-min(full_mdl,[],2));

%plot engagement versus PEV | motif marker by rec
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)
    idx = 1:14;    
    plot(engment_norm(j,:),full_mdl_norm(j,:),'marker',mark{j},'markersize',5,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (norm)')
ylabel('Performance (norm)')
P = polyfit(engment_norm(:),full_mdl_norm(:),1);
yfit = P(1)*engment_norm(:)+P(2);
plot(engment_norm(:),yfit,'r-','linewidth',2);
[rho,p] = corr(engment_norm(:),full_mdl_norm(:),'tail','right');
title(sprintf('More engagement = more predictability | Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')
xlim([-0.25 1.25]); ylim([-0.25 1.25])
%is this true for all brain regions? (yes)


%% Next we aimed to better understand this interactions. One hypothesis is that when a brain region is engaged with other areas it uses a subspace. 
%If so, then the majority of the variance should be captured by relatively
%few dimensions... To test this we used rRRR (show those plots). We found
%that, on average XX dimensions were needed to capture the vast majority of
%shared variances between an all other recorded regions. 
%For example, across recordings motif 6 needed XX dimensions
figure; hold on; 
ndim = 15;
thresh = 0.80;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
rrr_mdl = LoadVariable(data,'rel_performance',area_name);
rrr_mdl = squeeze(rrr_mdl(:,m,1:15));
arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(rrr_mdl,1));
plot(1:size(rrr_mdl,2), nanmean(rrr_mdl),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);  
d = arrayfun(@(n) find(rrr_mdl(n,:)>=thresh,1,'first'),1:size(rrr_mdl,1));
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel('performance');
title(sprintf('Motif %d needed %d to %d dim to capture %0.1d%% PEV',m,min(d), max(d),thresh*100),'Fontweight','normal')

%This was not due to our regression approach.
figure; hold on; 
cca_d = LoadVariable(datacca,'cca_dim',area_name);
rrr_d = LoadVariable(data,'rrr_dim',area_name);
%first plot our specific motif
plot([1,2],[rrr_d(:,m),cca_d(:,m)],'linewidth',0.5,'color',[0.75 0.75 0.75],'marker','o','markersize',4,'MarkerFaceColor',[0.75 0.75 0.75]);
errorbar(1,nanmean(rrr_d(:,m)),sem(rrr_d(:,m)),'linewidth',1.5,'color',[0.80,0.1,0.1],'marker','o','markersize',10);
errorbar(2,nanmean(cca_d(:,m)),sem(cca_d(:,m)),'linewidth',1.5,'color',[0.80,0.1,0.1],'marker','o','markersize',10);
[~,p] = ttest(rrr_d(:,m),cca_d(:,m),'tail','left');
set(gca,'xlim',[0.5 2.5],'ylim',[1 12],'xtick',[1,2],'xticklabel',{'rrr','cca'})
title(sprintf('Low dimensionality is not due to anlaysis method | Motif %d',m),'Fontweight','normal')
xlabel(sprintf('Method p=%0.2f',p));
ylabel('# of predictive dimensions');

% %Does dimensionality scale with activity
% engment_norm = (engment-min(engment,[],2))./(max(engment,[],2)-min(engment,[],2)); %normalize by the max prediction and activity across motifs per recording
% rrr_d_norm =  (rrr_d-min(rrr_d,[],2))./(max(rrr_d,[],2)-min(rrr_d,[],2));
% plot(engment_norm(:),rrr_d_norm(:),'o')
% close all
% [rho,p] = corr(engment(:),rrr_d_norm(:),'type','kendall')


%% Is dimensionality of subspace reciprical? 
%plot the reciprocal dimensionality
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










