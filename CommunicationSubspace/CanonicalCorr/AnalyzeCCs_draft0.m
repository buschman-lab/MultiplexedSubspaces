function AnalyzeCCs_draft0()
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

%remove noise motif and null motif
for i = 1:6
   idx = ismember(arrayfun(@(n) datacca{i}(n).motif, 1:size(datacca{i},2)),[2,16]);
   datacca{i}(idx)=[];   
   idx = ismember(arrayfun(@(n) data{i}(n).motif, 1:size(data{i},2)),[2,16]);
   data{i}(idx)=[];     
end


%% Question 1: When brain area is engaged by a motif, does it have more predictable varaibilty with other areas in the brian? 
m = 5; %[5,7,8];
close all
mout =11;
area_name = 'VIS';
% get performance of full model on visual area
full_mdl = [];
ex_mdl = []; %example recording
for i = 1:numel(data)
    for j = 1:size(data{1},2)
        idx = strcmp(data{i}(j).area_label,area_name);
        full_mdl(i,j) = cellfun(@(x) nanmean(bestLambda(x,1)), data{i}(j).cvl_ridge(idx),'UniformOutput',1);
        if i==1
            [~,temp] = cellfun(@(x) bestLambda(x,1), data{i}(j).cvl_ridge(idx),'UniformOutput',0);            
            ex_mdl(:,j) = 1-data{i}(j).cvl_ridge{idx}(:,temp{1});
        end
    end
end

% get engagement of visual areas across motifs and recordings
engment = NaN(numel(data),14);
for i = 1:numel(data)
    [yall, ~] = arrayfun(@(n) GetLocalTrialAverageStrength(data{i}(n).area_val,data{i}(n).area_label,'max'),1:size(data{i},2),'UniformOutput',0); 
    idx = strcmp(data{i}(1).area_label,area_name);
    engment(i,:) = cellfun(@(x) x(idx), yall);
end

%plot example recording Mean +/- SEM for full model 
figure; hold on; 
errorbar(1, nanmean(ex_mdl(:,m)),sem(ex_mdl(:,m)),'linestyle','none','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
plot(cat(1,ones(1,size(ex_mdl,1)),2*ones(1,size(ex_mdl,1))),ex_mdl(:,[m,mout])','marker','o','markersize',5,'color',[0.65 0.65 0.65],'linestyle','-')
errorbar(2,nanmean(ex_mdl(:,mout)),sem(ex_mdl(:,mout)),'linestyle','none','marker','o','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
set(gca,'xlim',[0.5 2.5],'ylim',[0.07 0.12],'xtick',[1,2],'xticklabel',[m,mout])
xlabel('motif'); ylabel('performance (10fold xval)');
title(sprintf('Example Performance motif %d area %s',m,area_name),'fontweight','normal')

%'vis' motifs should have more activity in it than, say 'non-visual' motif X
figure; hold on; 
bar(1,nanmean(engment(:,m),'all'),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0)
bar(2,nanmean(engment(:,mout),'all'),'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8],'EdgeAlpha',0)
for i = 1:size(engment,1)
   plot([1,2],[nanmean(engment(i,m)),nanmean(engment(i,mout))],'color',[0.25 0.25 0.25],'marker','o','markersize',5)
end
[~,p] = ttest(nanmean(engment(:,m),2),nanmean(engment(:,mout),2),'tail','right');
ylabel('Activity (normalized FR)')
xlabel(sprintf('motif p=%0.2f',p));
title(sprintf('Average "Engamgent" of %s area',area_name),'fontweight','normal')
set(gca,'xtick',1:2);


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
set(gca,'xtick',1:2);

%make engment relative
engment_z = zscore(engment,0,2);
full_mdl_z = zscore(full_mdl,0,2);

%plot engagement versus PEV | motif marker by rec
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)
    idx = 1:14;    
    plot(engment_z(j,:),full_mdl_z(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel('Performance (zscore)')
[rho,p] = corr(engment_z(:),full_mdl_z(:),'tail','right');
title(sprintf('More engagement = more predictability | Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')

%is this true for all brain regions? (yes)


%% Next we aimed to better understand this interactions. One hypothesis is that when a brain region is engaged with other areas it uses a subspace. 
%If so, then the majority of the variance should be captured by relatively
%few dimensions... To test this we used rRRR (show those plots). We found
%that, on average XX dimensions were needed to capture the vast majority of
%shared variances between an all other recorded regions. 

%plot the full example of the ridge and the RRR

%For example, across recordings motif 6 needed XX dimensions
figure; hold on; 
ndim = 15;
thresh = 0.80;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
ex_mdl=[];
d = NaN(1,numel(data));
for i = 1:numel(data)
    idx = strcmp(data{i}(j).area_label,area_name);
    if sum(idx)>0
        ex_mdl(i,:) = nanmean((1-data{i}(m).cvl_rrr{idx})/nanmean(bestLambda(data{i}(m).cvl_ridge{idx},0)));
        plot(1:size(ex_mdl,2), ex_mdl(i,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4);    
        d(i) = find(ex_mdl(i,:)>=thresh,1,'first');
    else
    end
end
plot(1:size(ex_mdl,2), nanmean(ex_mdl),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);    
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel('performance');
title(sprintf('Motif %d needed %d to %d dim to capture %0.1d%% PEV',m,min(d), max(d),thresh*100),'Fontweight','normal')

%This was not due to our regression approach... cca analysis...



%% Do the beta weights match what one might expect?
%plot the mean beta weight for each brain location for example motif
%another are in motif 5 = RSP ... so RSP and VIS should be closely correlated
%an example motif that only weakly engages VIS is motif 6 and 2
close all
m=8; %for something different use 2 or 11 and SS
area_name = 'VIS';
cur_d = 1; 
b = cell(1,numel(data));
for i = 1:numel(data)
    idx = strcmp(data{i}(1).area_label,area_name);
    if sum(idx)>0
        b{i}(:,1) = data{i}(m).rrr_B{idx}(:,cur_d); %just look at first dimension right now
        b{i}(:,2) = data{i}(m).grouping{idx};
    else
    end
end

%plot example recording
ugrp = cell(1,numel(data));
avg_beta = cell(1,numel(data));
contrib = cell(1,numel(data));
czscore = cell(1,numel(data));
for cur_rec = 1:numel(data)
    idx = strcmp(data{cur_rec}(1).area_label,area_name);
    if sum(idx)>0
        ugrp{cur_rec} = unique(b{cur_rec}(:,2));
        avg_beta = NaN(1,numel(ugrp));       
        [czscore{cur_rec},contrib{cur_rec}] = GetLocalTrialAverageStrength(data{cur_rec}(m).area_val,data{cur_rec}(m).area_label,'max',1);
    end       
end

figure; hold on; 
cur_rec=1;
for i = 1:numel(ugrp{cur_rec})
    tempidx = b{cur_rec}(:,2)==ugrp{cur_rec}(i);
    errorbar(ugrp{cur_rec}(i),nanmean(abs(b{cur_rec}(tempidx,1))),sem(abs(b{cur_rec}(tempidx,1))),'linestyle','none','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
%     errorbar(ugrp{cur_rec}(i),nanmean((b{cur_rec}(tempidx,1))),sem((b{cur_rec}(tempidx,1))),'linestyle','none','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
    bar(ugrp{cur_rec}(i),nanmean(abs(b{cur_rec}(tempidx,1))),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0);    
%     bar(ugrp{cur_rec}(i),nanmean((b{cur_rec}(tempidx,1))),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0);    
    avg_beta(i) = nanmean(abs(b{cur_rec}(tempidx,1)));
%     avg_beta(i) = nanmean((b{cur_rec}(tempidx,1)));
end
title(sprintf('Area %s motif %d ', area_name,m),'fontweight','normal');
ylabel('Beta weight (dimension 1)'); xlabel('brain area')
set(gca,'xtick',1:numel(data{cur_rec}(m).area_label),'xticklabel',data{cur_rec}(m).area_label);

%correlate with the relative engagement of that area for that recording
figure; hold on; 
temp = GetLocalTrialAverageStrength(data{cur_rec}(m).area_val,data{cur_rec}(m).area_label,'max',0);
plot(temp(1,ugrp{cur_rec}),avg_beta,'Marker','o','markersize',5,'linestyle','none')
x = temp(1,ugrp{cur_rec});
y = avg_beta;
P = polyfit(x(:),y(:),1);
yfit = P(1)*x(:)+P(2);
plot(x(:),yfit,'r-','linewidth',2); 
xlabel('Engagement of Visual Regions (not zscore)')
ylabel('Average |Beta|')
[rho,p] = corr(temp(1,ugrp{cur_rec})',avg_beta(:),'tail','right','type','pearson','rows','complete');
title(sprintf('Predictability vs Engagement Example Recording \n Motif %d | Rho=%0.2f p=%0.2f',m,rho,p),'fontweight','normal')

% confirm the above across all recordings. 
avg_beta = cell(1,numel(data));
avg_contrib = cell(1,numel(data));
for cur_rec = 1:numel(data)
%     temp = cellfun(@(x) nanmean(x), contrib{cur_rec});
    temp =czscore{cur_rec};
    btemp = zscore(b{cur_rec}(:,1));
%     btemp = b{cur_rec}(:,1);
    for i = 1:numel(ugrp{cur_rec})
        tempidx = b{cur_rec}(:,2)==ugrp{cur_rec}(i);
        %when comparing across recordings. beta weights should be zscored      
        avg_beta{cur_rec}(i) = nanmean(abs(btemp(tempidx))); 
%         avg_beta{cur_rec}(i) = nanmean((b{cur_rec}(tempidx,1))); 
        avg_contrib{cur_rec}(i) = temp(ugrp{cur_rec}(i));
    end   
end
y = cat(2,avg_beta{:});
x = cat(2,avg_contrib{:});
figure; hold on; 
plot(x,y,'Marker','o','markersize',5,'linestyle','none')
P = polyfit(x(:),y(:),1);
yfit = P(1)*x(:)+P(2);
plot(x(:),yfit,'r-','linewidth',2);    
xlabel('Engagement of Visual Regions (zscore)')
ylabel('Average |Beta| (zscore)')
[rho,p] = corr(x(:),y(:),'tail','right','type','pearson','rows','complete');
title(sprintf('Predictability vs Engagement All Recordings \n Motif %d | Rho=%0.2f p=%0.2f',m,rho,p),'fontweight','normal')

%% Engaging the same area in multiple ways = different beta weights. 
%across recordings, compare the engagement of motif X and Y across areas
m = 5;
mout=8;
area_name = 'VIS';
czscore = cell(1,numel(data));
b = cell(1,numel(data));
cur_d=1;
all_areas = data{1}(1).area_label;
for cur_rec = 1:numel(data)
    czscore{cur_rec} = NaN(2,numel(all_areas));
    b{cur_rec} = NaN(2,numel(all_areas)-1);
    idx = strcmp(data{cur_rec}(1).area_label,area_name);
    tempidx = ismember(all_areas,data{cur_rec}(1).area_label);
    if sum(idx)>0     
        czscore{cur_rec}(1,tempidx) = GetLocalTrialAverageStrength(data{cur_rec}(m).area_val,data{cur_rec}(m).area_label,'max',0);
        czscore{cur_rec}(2,tempidx) = GetLocalTrialAverageStrength(data{cur_rec}(mout).area_val,data{cur_rec}(mout).area_label,'max',0);
%         temp = zscore(data{cur_rec}(m).rrr_B{idx}(:,cur_d));
        temp = abs(zscore(data{cur_rec}(m).rrr_B{idx}(:,cur_d)));
        grp = data{cur_rec}(m).grouping{idx};
        b{cur_rec}(1,tempidx(1:end-1)) = arrayfun(@(n) nanmean(temp(grp==n)), unique(grp),'UniformOutput',1);
        temp = abs(zscore(data{cur_rec}(mout).rrr_B{idx}(:,cur_d)));
        grp = data{cur_rec}(mout).grouping{idx};
        b{cur_rec}(2,tempidx(1:end-1)) = arrayfun(@(n) nanmean(temp(grp==n)), unique(grp),'UniformOutput',1);
    end       
end
czscore = cat(1,czscore{:});
x = czscore(1:2:end,:);
y = czscore(2:2:end,:);
figure; hold on; 
errorbar(0.5:size(x,2)-0.5,nanmean(x),sem(x),'linestyle','none','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
errorbar(1:size(y,2),nanmean(y),sem(y),'linestyle','none','marker','none','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
bar(0.5:size(x,2)-0.5,nanmean(x),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0,'BarWidth',0.4);   
bar(1:size(y,2),nanmean(y),'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8],'EdgeAlpha',0,'BarWidth',0.4);   
for i = 1:size(x,2)
   plot([i-0.5,i],[x(:,i),y(:,i)],'marker','o','markersize',4,'color',[0.5 0.5 0.5])   
   if nanmean(x(:,i))>nanmean(y(:,i))
       [~,p] = ttest(x(:,i),y(:,i),'tail','right');
   else
       [~,p] = ttest(x(:,i),y(:,i),'tail','left');       
   end
   if p<0.05
      text(i-0.5,max(cat(1,x(:,i),y(:,i))),sprintf('%0.2f',p))
   end       
end

b = cat(1,b{:});
x = b(1:2:end,:);
y = b(2:2:end,:);
figure; hold on; 
errorbar(0.5:size(x,2)-0.5,nanmean(x),sem(x),'linestyle','none','marker','none','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
errorbar(1:size(y,2),nanmean(y),sem(y),'linestyle','none','marker','none','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
bar(0.5:size(x,2)-0.5,nanmean(x),'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1],'EdgeAlpha',0,'BarWidth',0.4);   
bar(1:size(y,2),nanmean(y),'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8],'EdgeAlpha',0,'BarWidth',0.4);   
for i = 1:size(x,2)
   plot([i-0.5,i],[x(:,i),y(:,i)],'marker','o','markersize',4,'color',[0.5 0.5 0.5])   
   if nanmean(x(:,i))>nanmean(y(:,i))
       [~,p] = ttest(x(:,i),y(:,i),'tail','right');
   else
       [~,p] = ttest(x(:,i),y(:,i),'tail','left');       
   end
   text(i-0.5,max(cat(1,x(:,i),y(:,i))),sprintf('%0.2f',p))       
end






%% compare the beta weights across dimensions



%


%% Are the beta weights being dominated by single units (plot across dimensions per source area and 1 target). 
%for motif that engages visual and area that's related, plot a PCA like
%plot of the relative contribution of each beta. 
close all
m=5;
area_name = 'VIS';
cur_d = 1; 
b = cell(1,numel(data));
for i = 1:numel(data)
    idx = strcmp(data{i}(1).area_label,area_name);
    if sum(idx)>0
        temp = data{i}(m).rrr_B{idx}(:,cur_d);
        grp = data{i}(m).grouping{idx};
        ugrp = unique(grp);
        temp = arrayfun(@(n) temp(grp==n),ugrp,'UniformOutput',0);
        temp = cellfun(@(x) sort(abs(x),'descend')/sum(abs(x)),temp,'UniformOutput',0);
        %get relative contribution
        b{i}(:,1) = cat(1,temp{:}); %just look at first dimension right now
        b{i}(:,2) = data{i}(m).grouping{idx};      
    else
    end
end

%plot the RSP for example recording
% figure; hold on; 
source_area = 'RSP';
cur_rec = 1; idx = strcmp(data{cur_rec}(1).area_label,source_area);
temp = b{cur_rec}(b{cur_rec}(:,2)==find(idx==1),1);
plot(cumsum(temp));

%across all regions, plot the fraction of the recorded neurons to capture 50% of
%the total beta weigths for that region.
close all
m=2;
area_name = 'VIS';
cur_d = 1; 
bfract = cell(1,numel(data));
for i = 1:numel(data)
    idx = strcmp(data{i}(1).area_label,area_name);
    if sum(idx)>0
        temp = data{i}(m).rrr_B{idx}(:,cur_d);
        grp = data{i}(m).grouping{idx};
        ugrp = unique(grp);
        temp = arrayfun(@(n) temp(grp==n),ugrp,'UniformOutput',0);
        temp = cellfun(@(x) sort(abs(x),'descend')/sum(abs(x)),temp,'UniformOutput',0);
        bfract{i} = cellfun(@(x) find(cumsum(x)>0.5,1,'first')/numel(x),temp,'UniformOutput',1);        
    else
    end
end
nanmean(cat(1,bfract{:}))
sem(cat(1,bfract{:}))

%% These results generalize across all areas and motifs
%see above for the example with another motif and two brain areas. 
%For example, do another 
figure; hold on; 
ndim = 15;
thresh = 0.80;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
ex_mdl=[];
d = NaN(1,numel(data));
for i = 1:numel(data)
    idx = strcmp(data{i}(j).area_label,area_name);
    if sum(idx)>0
        ex_mdl(i,:) = nanmean((1-data{i}(m).cvl_rrr{idx})/nanmean(bestLambda(data{i}(m).cvl_ridge{idx},0)));
        plot(1:size(ex_mdl,2), ex_mdl(i,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4);    
        d(i) = find(ex_mdl(i,:)>=thresh,1,'first');
    else
    end
end
plot(1:size(ex_mdl,2), nanmean(ex_mdl),'linestyle','-','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',5,'MarkerFaceColor',[0.8 0.1 0.1]);    
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel('performance');
title(sprintf('Motif %d needed %d to %d dim to capture %0.1d%% PEV',m,min(d), max(d),thresh*100),'Fontweight','normal')

%This was not due to our regression approach... cca analysis...



%% This suggests, a communication subspace... 
% to better undestand this; we focused on three motifs that engaged many of
% the same brain regions, but in different ways. 

%what I would expect is that there is significant differences in the
%incoming infromation but represented in either 1) the same neural
%population or 2) different neural populations. And that this similarity in
%the source region 
r = 1;
m = [5,7,8];
%compare the Vs
figure; hold on;
subplot(2,1,1); hold on; 
for i = m
plot(data{1}(i).rrr_V{r}(:,1));
end
corr(data{1}(5).rrr_V{r}(:,1),data{1}(7).rrr_V{r}(:,1))
m = [2,6,11];
%compare the Vs
subplot(2,1,2); hold on;
for i = m
plot(data{1}(i).rrr_V{r}(:,1));
end
corr(data{1}(5).rrr_V{r}(:,1),data{1}(11).rrr_V{r}(:,1))
corr(data{1}(6).rrr_V{r}(:,1),data{1}(2).rrr_V{r}(:,1))





%% first, we compared their subspace dimensionalities but in agnostic to an arbitraty cutoff like above. 
%first we asked if they have different subspaces dimensionality within visual?
%this can be done as a cross validation
m = [5,7,8];
temp = nchoosek(m,2);
vis=[];
for i = 1:size(temp,1)
    [vis,pval,rr] = SubspaceROC_single(data,'VIS',temp(i,:),1);
end
%i.e. no different when they engage the same region
%what what about when they engage different regions? 
mos=[];
for i = 1:size(temp,1)
    [mos,pval,rr] = SubspaceROC_single(data,'MOs',temp(i,:),1);
end

% motifs that engaged a brain region more exhibit lower dimensionality than other motifs 
temp = combvec(1:14,1:14)';
auc=NaN(14,14,numel(data));
for i = 1:size(temp,1)
    if temp(i,1)-temp(i,2)==0
        auc(temp(i,1),temp(i,2),:)=NaN(1,numel(data));
    else
        auc_temp = SubspaceROC_single(data,area_name,temp(i,:),0);
        auc(temp(i,1),temp(i,2),1:size(auc_temp,2))=auc_temp;
    end
end
%get the column
temp = squeeze(nanmean(auc,1))';
temp = zscore(temp,[],2);
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)   
    plot(engment_z(j,:),temp(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel(sprintf('# of dimensions to explain %0.2d%% PEV',thresh*100))
[rho,p] = corr(engment_z(:),temp(:),'tail','left','type','pearson','rows','complete');
title(sprintf('Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')
%this is still true when controlling for the 
% [rho,pval] = partialcorr(engment_z(:),temp(:),full_mdl_z(:),'Type','Pearson','tail','left');

%% does this finding generalize (No) ... may want to reverse order with the subsequent finding. 
area_name = 'SSp-bfd';
% get performance of full model on visual area
full_mdl = NaN(numel(data),size(data{1},2));
for i = 1:numel(data)
    for j = 1:size(data{1},2)
        idx = strcmp(data{i}(j).area_label,area_name);
        if sum(idx)>0
            full_mdl(i,j) = cellfun(@(x) nanmean(bestLambda(x,1)), data{i}(j).cvl_ridge(idx),'UniformOutput',1);
        end
    end
end

% get engagement of visual areas across motifs and recordings
engment = NaN(numel(data),14);
for i = 1:numel(data)
    [yall, ~] = arrayfun(@(n) GetLocalTrialAverageStrength(data{i}(n).area_val,data{i}(n).area_label,'max'),1:size(data{i},2),'UniformOutput',0); 
    idx = strcmp(data{i}(1).area_label,area_name);
    if sum(idx)>0
        engment(i,:) = cellfun(@(x) x(idx), yall);
    end
end

full_mdl = zscore(full_mdl,[],2);
engment = zscore(engment,[],2);

% motifs that engaged a brain region more exhibit lower dimensionality than other motifs 
temp = combvec(1:14,1:14)';
auc=NaN(14,14,numel(data));
for i = 1:size(temp,1)
    if temp(i,1)-temp(i,2)==0
        auc(temp(i,1),temp(i,2),:)=NaN(1,numel(data));
    else
        auc(temp(i,1),temp(i,2),:) = SubspaceROC_single(data,area_name,temp(i,:),0);
    end
end
%get the column
temp = squeeze(nanmean(auc,1))';
temp = zscore(temp,[],2);
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)   
    plot(engment_z(j,:),temp(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel(sprintf('# of dimensions to explain %0.2d%% PEV',thresh*100))
[rho,p] = corr(engment_z(:),temp(:),'tail','left','type','pearson','rows','complete');
if rho>0
    [rho,p] = corr(engment_z(:),temp(:),'tail','right','type','pearson','rows','complete');
else
    [rho,p] = corr(engment_z(:),temp(:),'tail','left','type','pearson','rows','complete');
end


title(sprintf('Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')
%this is still true when controlling for the 
% [rho,pval] = partialcorr(engment_z(:),temp(:),full_mdl_z(:),'Type','Pearson','tail','left','Rows','complete')



%% This effectively argues for this idea of communication subspaces



% However, this dependent on the brain regions



%in other words, when vis is strongly engaged and interacting with other 
%areas, it's local dimensionality ___ and it's
%its communication dimensionality decreases.... i.e. supporting this idea
%that 



%Next we were interested in how this compared to other motifs that shared
%the same area. 
[auc,pval,rr] = SubspaceROC_single(data,area_name,[m,mout],1);



%in other words, when vis is strongly engaged and interacting with other 
%areas, it's local dimensionality ___ and it's
%its communication dimensionality decreases.... i.e. supporting this idea
%that 

%Next we were interested in how 


%Interestingly the dimensionality of this differed across motifs. To do this we used our ROC.. This
%limits arbitrary selection of the number of dimensions. 

% Interestingly, when adjusted for the higher amount of variance explained, area engaged more explained more variance in less dimensions.. consistent with this idea of a communication subspace. 
%for example: show a pairwise comparison. 

%%on potential interesting points that would support this idea of subspaces
%is if the partial correlation between 
d = NaN(numel(data),size(data{1},2));
for i = 1:numel(data)
    for j = 1:size(data{i},2)
        idx = contains(data{i}(j).area_label,area_name);
        d(i,j) = find(nanmean((1-data{i}(j).cvl_rrr{idx})/nanmean(bestLambda(data{i}(j).cvl_ridge{idx},0)))>=thresh,1,'first');
    end
end
d = zscore(d,0,2);
[rho,pval] = partialcorr(engment_z(:),d(:),full_mdl_z(:),'Type','Pearson','tail','left');



%% We were interested in whether these findings generalized across brain areas


%% When a brain region is engaged by multiple motifs 
%plot example 
figure; hold on; 
cur_rec=1
idx = contains(data{i}(j).area_label,area_name);
ex_mdl = (1-data{cur_rec}(m).cvl_rrr{idx})/nanmean(bestLambda(data{cur_rec}(m).cvl_ridge{idx},1));
errorbar(1:size(ex_mdl,2), nanmean(ex_mdl),sem(ex_mdl),'linestyle','none','marker','o','color',[0.8 0.1 0.1],'linewidth',1.5,'markersize',10);
ex_mdl = (1-data{cur_rec}(mout).cvl_rrr{idx})/nanmean(bestLambda(data{cur_rec}(mout).cvl_ridge{idx},1));
errorbar(1:size(ex_mdl,2), nanmean(ex_mdl),sem(ex_mdl),'linestyle','none','marker','o','color',[0.1 0.1 0.8],'linewidth',1.5,'markersize',10);
xlim([0.5 10]);

%compare with ROC
[auc,pval,rr] = SubspaceROC_single(data,area_name,[m,mout],1);

%Across recordings | number of dimensions needed to reach x threshold
a = squeeze(rr(:,1,:));
b = squeeze(rr(:,2,:));
figure; hold on; 
plot(a','linewidth',0.5,'color',[0.8 0.1 0.1])
plot(b','linewidth',0.5,'color',[0.1 0.1 0.8])

%plot engagement versus dimensionality (to reach XX percent) 
full_mdl = [];
for i = 1:numel(data)
    for j = 1:size(data{1},2)
        idx = contains(data{i}(j).area_label,area_name);
        full_mdl(i,j) = find(nanmean((1-data{i}(j).cvl_rrr{idx})/nanmean(bestLambda(data{i}(j).cvl_ridge{idx},1)))>=0.75,1,'first');
    end
end

mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark) 
    plot(engment_z(j,:),full_mdl(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel('Number of dimensions')
[rho,p] = corr(engment_z(:),full_mdl(:),'tail','right','type','kendall');
title(sprintf('Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')

%repeat using number of CCA dimensions



%plot engagement versus ROC across all motifs
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)
    idx = 1:14;    
    plot(engment_z(j,:),full_mdl_z(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel('Performance (zscore)')
[rho,p] = corr(engment_z(:),full_mdl_z(:),'tail','right');
title(sprintf('Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')

%other questions: could also do relative to the dimensionality of that
%region during that motif. I.e. is lower D just because the dimensionality
%is lower D? or is it actually the opposite? 

%% Furthermore, these subspaces are different (from beta weights). 


% do more similarly engagement motifs have more similar subspaces? 



%% Choose 3 Motifs that engage the same brain areas in different ways
m = [5,7,8]; %all engage rsp and vis regions in different ways
region = [4,8];

%INCLUDE A REGION WITH LITTLE TO NO VISUAL CORTICAL ACTIVITY

% Are they using subspaces? 
%general subspace plot
%plot of the local dimensionality


% Do they differ in their subspace dimensionality? 

[auc_all,pval_all] = SubspaceROC(data_rrr(m)); close all;
auc = cellfun(@(x) nanmean(x,3), auc_all(region),'UniformOutput',0);
auc = cat(3,auc{:});
for cur_a = region
    data = arrayfun(@(n) (1-data_rrr(n).cvl_rrr{cur_a})/max((1-data_rrr(n).cvl_rrr{cur_a}),[],'all'), m,'UniformOutput',0);
    %plot example curve
    figure; hold on; 
    t = tiledlayout(1,nchoosek(numel(m),2));
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    count = 1;
    for i = 1:size(data,2)
        for j = 1:size(data,2)               
            if j>i
                nexttile
                hold on;
                arrayfun(@(n) plot([0,data{i}(n,:),1],[0,data{j}(n,:),1],'linewidth',0.25,'color',[0.75 0.75 0.75]),1:10)
                plot([0,nanmean(data{i}),1],[0,nanmean(data{j}),1],'linestyle','-','marker','none','linewidth',2,'color','k'); hold on;        
                xlabel(sprintf('Motif %d',i));
                ylabel(sprintf('Motif %d',j));
                plot(0:1,0:1,'linewidth',1,'color','r'); 
                axis square
                count = count+1;
                title(sprintf('auc=%0.2f',nanmean(auc_all{cur_a}(i,j,:))),'fontweight','bold')
            end        
        end 
    end
    title(t,area_label{cur_a})
    set(gcf,'units','normalized','position',[0 0 1 1])
end
%answer yes.... but it's region specific -- could combine across recs...
%potentially

%how does this relate to how strongly they engage a region? 
[yall, ~] = arrayfun(@(n) GetLocalTrialAverageStrength(data_rrr(n).area_val,data_rrr(n).area_label,'max'),m,'UniformOutput',0);
yall = cat(1,yall{:});
yall = yall(:,region);

%for example, do motifs that engage a brain area more similarly, have similar dimensionality? 
ydiff = NaN(size(yall,1),size(yall,1),numel(region));
for i = 1:size(ydiff,1)
   for j = 1:size(ydiff,2)
       if i~=j
           ydiff(i,j,:) = yall(i,:)-yall(j,:);
       else
           ydiff(i,j,:) = NaN(1,size(ydiff,3));
       end
   end
end

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(region)
    nexttile; hold on;
    a = triu(ydiff(:,:,i),1); %only use the upper since symmetric
    a(a==0)=[];
    b = triu(auc(:,:,i),1); 
    b(b==0)=[];
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('\Delta Activity (neg=M1<active than M2)');
    ylabel('\Delta ROC (low: M1=<dim than M2)');
    title(['Area ',area_label{region(i)} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end

% In regions where they have the same subspace dimensionality... are they
% using the same subspace? How does this compare to regions where they have different subspace
% dimensionality? 

[rsq_full,rsq_byarea] = BetaGeneralization(data_rrr(m),0);
rsq_full = cat(3,rsq_full{region});
rsq_byarea = cat(1,rsq_byarea{region});

%is the interaction between the two brain regions dominated by the same
%region 1 to 2 and 2 to 1
rsq_byarea = rsq_byarea(:,region);

%if yes = shared representations
%if no = unique subspaces











%%
% motifs that engaged a brain region more exhibit lower dimensionality than other motifs 
temp = combvec(1:14,1:14)';
auc=NaN(14,14,numel(data));
for i = 1:size(temp,1)
    if temp(i,1)-temp(i,2)==0
        auc(temp(i,1),temp(i,2),:)=NaN(1,numel(data));
    else
        auc(temp(i,1),temp(i,2),:) = SubspaceROC_single(data,area_name,temp(i,:),0);
    end
end
%get the column
temp = squeeze(nanmean(auc,1))';
% temp = zscore(temp,[],2);
mark = {'+','o','*','x','v'};
figure; hold on; 
for j = 1:numel(mark)   
    plot(engment_z(j,:),temp(j,:),'marker',mark{j},'markersize',3,'color','k','linestyle','none')
end
xlabel('Engagement of Visual Regions (zscore)')
ylabel(sprintf('# of dimensions to explain %0.2d%% PEV',thresh*100))
[rho,p] = corr(engment_z(:),temp(:),'tail','left','type','pearson');
title(sprintf('Rho=%0.2f p=%0.2f',rho,p),'fontweight','normal')

















%% Plot the subspace ROC between motifs
[auc_all,pval_all] = SubspaceROC(data_rrr); close all;

auc = cellfun(@(x) nanmean(x,3),auc_all,'UniformOutput',0);
auc = cat(3,auc{:});

%% Get the activity in each brain region during a given motif
area_label = data_rrr(1).area_label;
[yall, ~] = arrayfun(@(n) GetLocalTrialAverageStrength(data_rrr(n).area_val,data_rrr(n).area_label,'max'),1:size(data_rrr,2),'UniformOutput',0);

yall = cat(1,yall{:});

%compare differences in engagement of a region to differences in the
%subspace structure; 


ydiff = NaN(size(yall,1),size(yall,1),numel(area_label));
for i = 1:size(ydiff,1)
   for j = 1:size(ydiff,2)
       if i~=j
           ydiff(i,j,:) = yall(i,:)-yall(j,:);
       else
           ydiff(i,j,:) = NaN(1,size(ydiff,3));
       end
   end
end

%plot matrix per brain area
%visualize the mean
for cur_a = 1:numel(area_label)
figure; hold on;
a = imagesc(ydiff(:,:,cur_a)); 
c=colorbar; colormap magma
ylabel(c,'/Delta norm FR | light: x motif higher activity than y motif')
set(gca,'xlim',[0.5 size(data_rrr,2)+0.5],'ylim',[0.5 size(data_rrr,2)+0.5]);
title(['Area ',data_rrr(1).area_label{cur_a}],'FontWeight','normal');
xlabel('"x" motif'); ylabel('"y" motif');
end


figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = triu(ydiff(:,:,i),1); %only use the upper since symmetric
    a(a==0)=[];
    b = triu(auc(:,:,i),1); 
    b(b==0)=[];
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('\Delta Activity (neg=M1<active than M2)');
    ylabel('\Delta ROC (low: M1=<dim than M2)');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end

%% to clarify, let's take one example motif in one brain region and compare to the other motifs in that brain region 

%brain area 3 motif 5 (vis) versus all other motifs 
m = 1;
cur_a = 3;

figure('position',[681         559        1174         420]); hold on; 
tiledlayout(1,3);
x = ydiff(m,:,cur_a);
nexttile; hold on;
bar(x); 
ylabel(['\Delta peak normalized FR',sprintf('(Motif%d - Motif X)',m)]); 
xlabel('motifs')
title('example differences in neural activity','FontWeight','normal')

y = auc(m,:,cur_a);
nexttile; hold on; bar(y); 
ylabel(['\Delta mean AUC',sprintf('(>0.5 = Motif%d lower dimensional)',m)]); 
xlabel('motifs')
title('example differences in dimensionality','FontWeight','normal')
plot([0 15],[0.5 0.5],'linewidth',2)
ylim([0.3 0.7])

%their correlation
nexttile; hold on; 
plot(x,y,'linestyle','none','marker','o','markersize',5);
xlabel('\Delta Activity');
ylabel('\Delta ROC');
x(isnan(x))=[];
y(isnan(y))=[];
rho = corr(x',y','rows','complete');
if rho<0
    [rho,pval] = corr(x',y','rows','complete','tail','left');
else
    [rho,pval] = corr(x',y','rows','complete','tail','right');
end
title(['Area ',area_label{cur_a} sprintf('  | Motif %d r=%0.2f p=%0.2f',m,rho,pval)],'FontWeight','normal');
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
plot(x,yfit,'r-','linewidth',2);    
yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
set(gca,'xlim',xvals,'ylim',yvals)
axis square
grid on


%% Question 1: do motifs that engage a brain region similarly (abs(delta activity)) exhibit similar dimensionality abs(auc))


figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = abs(triu(ydiff(:,:,i),1)); %only use the upper since symmetric
    a(a==0)=[];    
    b = abs(triu(auc(:,:,i)-0.5,1)); 
    b(b==0)=[];
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('simility in "engagement" of brain region');
    ylabel('Similarity in dimensionality');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end


%% Question 2: do motifs that more strongly engage a brain region do so with larger or smaller dimensionality? 

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = yall(:,i);    
    %average ROC for that motif relative to all other motifs (column so
    %lower = lower dimensional
    b = nanmean(auc(:,:,i),1)'; 
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('simility in "engagement" of brain region');
    ylabel('Similarity in dimensionality');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end

close all;

%% Question 2 in a different way: 
%do motifs that more strongly engage a brain region do so with larger or smaller dimensionality
% number of significant dimensions from CCA
r = cat(1,data_cca.r);
r = cellfun(@(x) numel(x), r,'UniformOutput',1);

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = yall(:,i);    
    %average ROC for that motif relative to all other motifs (column so
    %lower = lower dimensional
    b = r(:,i); 
    rho = corr(a(:),b(:),'rows','complete','type','spearman');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left','type','spearman');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right','type','spearman');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('Strength of "engagement" of brain region');
    ylabel('number of significant CCA dimensions');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end

%% compare activity to the local dimensionality
[~,~,d,~] = localDimensionality(data_fa);
d = d';

figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = yall(:,i);
    b = d(:,i);
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('Activity');
    ylabel('Dimensionality');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end



%% Do motifs with similar engagement of a brain region show similar representations or less similar? 

[rho_all,midx_all] = SubspaceCorr(data_rrr); 

% loop through subspace dimensions
% loop through brain areas
rho_temp = cellfun(@(x) x(:,:,2),rho_all,'UniformOutput',0);
rho_temp = cat(3,rho_temp{:});
rho_temp = fisherZ(rho_temp);


figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = abs(triu(ydiff(:,:,i),1)); %only use the upper since symmetric
    a(a==0)=[];
    b = triu(rho_temp(:,:,i),1); 
    b(b==0)=[];
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('\Delta Activity (neg=M1<active than M2)');
    ylabel('Similarity in beta weights');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end

title(t,'The larger the difference between how motifs engage a region, the less similarity in beta weights','VerticalAlignment','bottom')

%% now test if motifs that more strongly engage a brain region have more dissimilarity in their subspace relative to others
% loop through subspace dimensions
% loop through brain areas
for cur_d = 1:6
rho_temp = cellfun(@(x) x(:,:,cur_d),rho_all,'UniformOutput',0);
rho_temp = cat(3,rho_temp{:});
rho_temp = fisherZ(rho_temp);


figure('units','normalized','position',[0 0 1 1]); hold on; 
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(area_label)
    nexttile; hold on;
    a = yall(:,i);
    a(isnan(a))=[];
    %average similarity between dimension of motif and any other dimensions of other motifs (i.e. column)
    b = nanmean(rho_temp(:,:,i),1)';        
    rho = corr(a(:),b(:),'rows','complete');
    if rho<0
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','left');
    else
        [rho,pval] = corr(a(:),b(:),'rows','complete','tail','right');
    end
    plot(a(:),b(:),'linestyle','none','marker','o','markersize',3);
    xlabel('motif "engagement" of brain region');
    ylabel('Similarity in beta weights');
    title(['Area ',area_label{i} sprintf('  |  r=%0.2f p=%0.2f',rho,pval)],'FontWeight','normal');
    P = polyfit(a(:),b(:),1);
    yfit = P(1)*a(:)+P(2);
    plot(a(:),yfit,'r-','linewidth',2);    
    yvals = get(gca,'ylim'); xvals = get(gca,'xlim');
    plot([0,0],yvals,'linestyle',':','linewidth',1,'color','k')
    plot(xvals,[0.5 0.5],'linestyle',':','linewidth',1,'color','k')
    set(gca,'xlim',xvals,'ylim',yvals)
    axis square
    grid on
end
title(t,sprintf('Dimension %d',cur_d));

end

%%This analysis is what gets us the topography question. If we see that
%%in some brain area more similar motifs engage similar patterns
%%(negative correlation above) and stronger motifs engage the more similar
%%patterns, this suggests a... this gets us to whether they are using
%%'significantly similar subspaces'... while other areas use unique
%%subspaces 
%i.e. in some nueral populations there is significant mixing of
%representation or lack there off. Particularly Make a plot that is
%basically brain area by 'mixing' level (correlation in betas) across
%dimensions... maybe a 3D plot with the PEV of that dimension too... so
%those in which the subspace is dominanted by low dimension there is less
%mixing. Other regions there is more mixing. 
%... need to come up with a significance test for whether something is more
%mixed than by chance. 


%the final question in this across motifs analysis 
%is whether the motifs with similar beta weights really are using the same
%subspace... if so then project the activity of one should be able to
%explain the vast majority of the variance in the other. 





%you could then do the same analysis but for when brain area A communicates
%with brain area B or brain area C (but would require synethtic data). 














%% Compare local dimensionality and communication subspace dimensions and CCA dimensions
cca_d = cell(1,size(data_fa,2));
local_d = cell(1,size(data_fa,2));
pred_d = cell(1,size(data_fa,2));
pred_d_lower = cell(1,size(data_fa,2));
score = cell(1,size(data_fa,2));
pneu = cell(1,size(data_fa,2));
tneu = cell(1,size(data_fa,2));
e_pt = cell(1,size(data_fa,2));
e_pev = cell(1,size(data_fa,2));
for cur_fit = 1:size(data_fa,2)
    %number of significant CCA dimensions
    cca_d{cur_fit} = cellfun(@(x) numel(x),data_cca(cur_fit).r,'UniformOutput',1); 
    
    %number of dimensions to reach 90% of local variance
    idx = arrayfun(@(n) find(data_fa(cur_fit).paired_areas(:,2)==n,1,'first'),unique(data_fa(cur_fit).paired_areas(:,2)),'UniformOutput',1); 
    local_d{cur_fit} = cellfun(@(x) find((1-x)>=0.99,1,'first'),data_fa(cur_fit).cvl_fa_target(idx),'UniformOutput',0); 
       
    %# of dimensions needed to reach 90% of predictable variance in activty    
    score{cur_fit} = cellfun(@(x) nanmean(bestLambda(x)), data_rrr(cur_fit).cvl_ridge,'UniformOutput',0);
    [pneu{cur_fit},tneu{cur_fit}] = cellfun(@(x) size(x),data_rrr(cur_fit).rrr_b,'UniformOutput',0);
    pred_d{cur_fit} = cellfun(@(x,y) find(nanmean(1-x)>=(0.8*y), 1,'first'),data_rrr(cur_fit).cvl_rrr,score{cur_fit},'UniformOutput',0)';
    pred_d_lower{cur_fit} = cellfun(@(x,y) find(nanmean(1-x)>=(0.5*y), 1,'first'),data_rrr(cur_fit).cvl_rrr,score{cur_fit},'UniformOutput',0)';
    
    %get the elbow pt... lower elbow with higher captured variance means more subspaces
    e_pt{cur_fit} = cellfun(@(x) GetElbow(x),data_rrr(cur_fit).cvl_rrr,'UniformOutput',0); 
    e_pev{cur_fit} = cellfun(@(x,y) (1-nanmean(x(:,GetElbow(x))))./y,data_rrr(cur_fit).cvl_rrr,score{cur_fit},'UniformOutput',0);
    
    [pneu{cur_fit},tneu{cur_fit}] = cellfun(@(x) size(x),data_rrr(cur_fit).rrr_b,'UniformOutput',0);
    
end

pred_d = cat(1,pred_d{:});
pred_d_lower = cat(1,pred_d_lower{:});
%set empties to NaN
for i = 1:numel(pred_d); if isempty(pred_d{i}); pred_d{i} = NaN; end; end
for i = 1:numel(pred_d_lower); if isempty(pred_d_lower{i}); pred_d_lower{i} = NaN; end; end
pred_d = cell2mat(pred_d);
pred_d_lower = cell2mat(pred_d_lower);
pneu = cell2mat(cat(2,pneu{:}));
tneu = cell2mat(cat(2,tneu{:}));
local_d = cell2mat(cat(1,local_d{:}));
score = cell2mat(cat(2,score{:}))';
e_pev = cell2mat(cat(2,e_pev{:}));
e_pt = cell2mat(cat(2,e_pt{:}));
m = repmat({'+','o','*','x','v','d','^','s'},1,size(data_fa,2));
cca_d = cat(2,cca_d{:})';

figure; hold on; 
plot(local_d,pred_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
plot(1:max(cat(1,local_d,pred_d)),1:max(cat(1,local_d,pred_d)),'linestyle','--','linewidth',2,'color','k')
xlabel('# of local dimensions (99% of optimized... note pca is MUCH larger)') 
ylabel('# of pred dim for 75% of explainABLE variance') 
title('predictive versus local dim')

figure; hold on; 
histogram(local_d./pred_d,'BinWidth',0.1,'EdgeAlpha',0,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5)

%histogram of PEV explained by full model
figure; hold on; 
histogram(score,'BinWidth',0.01,'EdgeAlpha',0,'FaceColor',[0.8 0.1 0.1],'FaceAlpha',0.5)
xlabel('Performance (rsq) of the full model');
ylabel('# of Fits')
title('Amount of trial-by-trial variance explained by Ridge Regression Full Model','Fontweight','normal')

%impact of neuron counts
figure; hold on; 
t=tiledlayout(2,2);
nexttile
plot(pneu,local_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('local d'); xlabel('# pred neu');
title(sprintf('rho %0.2f',corr(pneu',local_d)))
nexttile
plot(tneu,local_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('local d'); xlabel('# targ neu');
title(sprintf('rho %0.2f',corr(pneu',local_d)))
nexttile
plot(pneu,score,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('Full Model performance'); xlabel('# pred neu');
title(sprintf('rho %0.2f',corr(pneu',score)))
nexttile
plot(tneu,score,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('Full Model  performance'); xlabel('# targ neu');
title(sprintf('rho %0.2f',corr(pneu',score)))
title(t,'Number of neurons matters')

%impact of neuron counts
figure; hold on; 
t=tiledlayout(2,2);
nexttile
plot(pneu,pred_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('pred d'); xlabel('# pred neu');
title(sprintf('rho %0.2f',corr(pneu',pred_d)))
nexttile
plot(tneu,pred_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('pred d'); xlabel('# targ neu');
title(sprintf('rho %0.2f',corr(tneu',pred_d)))
nexttile
plot(pneu,cca_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('CCA dimensions'); xlabel('# pred neu');
title(sprintf('rho %0.2f',corr(pneu',cca_d)))
nexttile
plot(tneu,cca_d,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
ylabel('CCA dimensions'); xlabel('# targ neu');
title(sprintf('rho %0.2f',corr(tneu',cca_d)))
title(t,'Number of neurons matters')


%relationship between elbow point and pt
figure; hold on; 
arrayfun(@(n) plot(e_pt(n),e_pev(n),'linestyle','none','marker',m{n},'markersize',5,'color',[0.8 0.1 0.1]),1:numel(m))
xlabel('# Dimensions at Elbow') 
ylabel('Fraction of Explainable Variance') 

%combined together
figure; hold on; 
arrayfun(@(n) plot(local_d(n)./pred_d(n),pred_d(n)./pred_d_lower(n),'linestyle','none','marker',m{n},'markersize',5,'color',[0.8 0.1 0.1]),1:numel(m))
xlabel('# local/pred') 
ylabel('# D for 90% / #D for 50%') 
title('Comparing Subspace D w/ "slope"')

%histogram of the number of dimensions at the elbow 
figure; hold on; histogram(e_pt(:))

%histogram of PEV at the elbow 
figure; hold on; histogram(e_pev(:))

figure; hold on; 
plot(e_pt,e_pev,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
xlabel('# of dimensions at elbox') 
ylabel('% explainable variance')


%comparing CCA and rRRR dimensions
figure; hold on; 
rng('default')
temp = rand(numel(cca_d),1)/2-0.25;
temp2 = rand(numel(cca_d),1)/2-0.25;
plot(cca_d+temp,pred_d+temp2,'linestyle','none','marker','o','markersize',5,'color',[0.8 0.1 0.1])
plot(1:max(cat(1,cca_d,pred_d)),1:max(cat(1,cca_d,pred_d)),'linestyle','--','linewidth',2,'color','k')
xlim([0 10]); ylim([0 10]);
ylabel('RRR Dimensions'); xlabel('CCA Dimensions')
[rho,p] = corr(cca_d,pred_d,'type','Pearson','tail','right');
title(sprintf('80%% Rho %0.2f, pval %0.2g',rho,p),'fontweight','normal')
P = polyfit(cca_d,pred_d,1);
yfit = P(1)*cca_d+P(2);
plot(cca_d,yfit,'r-.','linewidth',2);


%% Are subspaces different between motifs. 
%grab statistical comparison
s_all = cell(1,size(data_rrr,2));
for cur_fit = 1:size(data_rrr,2)
    %for a given target region
    score = cellfun(@(x) nanmean(bestLambda(x)), data_rrr(cur_fit).cvl_ridge,'UniformOutput',0);
    s = cellfun(@(x,y) (1-reshape(x',[1,size(x,2),size(x,1)]))./y, data_rrr(cur_fit).cvl_rrr_unique,score,'UniformOutput',0);
    s = cat(1,s{:});    
    s_all{cur_fit} = reshape(s,size(s,1)*size(s,2),size(s,3));
end
s_all = cat(3,s_all{:});

%get the anova
temp = NaN(size(s_all));
p = NaN(size(s_all,1),1);
for i = 1:size(s_all,1)
%     p(i) = anova1(squeeze(s_all(i,:,:)),[],'on');
    a = squeeze(s_all(i,:,:));       
    a = zscore(a(:));
    a = reshape(a,size(s_all,2),size(s_all,3));
    p(i) = anova1(a,[],'off');  
    temp(i,:,:) = a;
end

figure; hold on; 
imagesc(squeeze(nanmean(temp,2)));
ylabel('contribution of different pairs of regions')
xlabel('motif');
colorbar
title({'Subspace strength','Zscored (including cross validations'},'fontweight','normal'); 


%% Compare the number of dimensions needed 
pred_d = cell(1,size(data_rrr,2));
for cur_fit = 1:size(data_rrr,2)
    %for a given target region
    score = cellfun(@(x) nanmean(bestLambda(x)), data_rrr(cur_fit).cvl_ridge,'UniformOutput',0);
    pred_d{cur_fit} = cellfun(@(x,y) find(nanmean(1-x)>=(0.75*y), 1,'first'),data_rrr(cur_fit).cvl_rrr,score,'UniformOutput',0)';
end

pred_d = cat(2,pred_d{:});
%set empties to NaN
for i = 1:numel(pred_d); if isempty(pred_d{i}); pred_d{i} = NaN; end; end
pred_d = cell2mat(pred_d);

figure; hold on; 
imagesc(pred_d,[0 10]); colormap magma
title('number of dimensions per target area'); colorbar
ylabel('target area'); 
set(gca,'ytick',[1:8],'YTickLabel',data_rrr(1).area_label,'ylim',[0.5 8.5],'xlim',[0.5 14.5])
xlabel('motif');

%normalize to zero and on per row
pred_d = (pred_d-min(pred_d,[],2))./(max(pred_d,[],2)-min(pred_d,[],2));
figure; hold on; 
imagesc(pred_d); colormap magma
title('number of dimensions (normalized) per target area')
ylabel('target area'); 
xlabel('motif');
colorbar

%% Are the subspaces expected? 
%plot the variance in different see regions versus subspace strength


%% How are the subspaces represented? 
%get the mixing of beta weights

inweight = cell(1,size(data_rrr,2));
for cur_fit = 1:size(data_rrr,2)
    %rRRR keeps all beta weights - the rank limiting is on the V    
    score = cellfun(@(x) nanmean(bestLambda(x)), data_rrr(cur_fit).cvl_ridge,'UniformOutput',0);
    d = cellfun(@(x,y) find(nanmean(1-x)>=(0.75*y), 1,'first'),data_rrr(cur_fit).cvl_rrr,score,'UniformOutput',1)';
    B = arrayfun(@(n) data_rrr(cur_fit).rrr_B{n}(:,1:d(n)), 1:numel(d),'UniformOutput',0);
    grp = arrayfun(@(n) data_rrr(cur_fit).grouping{n},1:numel(d),'UniformOutput',0);
    [d,contrib,p] = cellfun(@(x,y) BetaMixing(x,y),B,grp,'UniformOutput',0);
    inweight{cur_fit} = cat(1,d{:});
end
inweight = cat(3,inweight{:});

figure; hold on; 
t = tiledlayout(1,3);
temp = {'% fully mixed (no sig anova)','% partial mixed (mutlcompare)','% solo contribution'};
for i = 1:3    
    nexttile
    imagesc(squeeze(inweight(:,i,:)))
    ylabel('area');
    set(gca,'ytick',1:8,'YTickLabel',data_rrr(1).area_label)
    xlabel('motif');
    title(gca,temp{i});
    colorbar
end
title(t,{'Mixing of representations within Beta Weights to a Target regions'}); 

% compare mixing per motif
% compare mixing across brain regions

%plot with predictive dimension
inweight = cell(1,size(data_rrr,2));
for cur_fit = 1:size(data_rrr,2)
    %rRRR keeps all beta weights - the rank limiting is on the V    
    temp = [];
    for j = 1:10
        B = arrayfun(@(n) data_rrr(cur_fit).rrr_B{n}(:,j), 1:8,'UniformOutput',0);
        grp = arrayfun(@(n) data_rrr(cur_fit).grouping{n},1:numel(d),'UniformOutput',0);
        [a,~,~] = cellfun(@(x,y) BetaMixing(x,y),B,grp,'UniformOutput',0);
        a = cat(1,a{:});
        temp(j,:) = nanmean(a,1);
    end
    inweight{cur_fit} = temp;
end
inweight = cat(3,inweight{:});

figure; hold on; 
shadedErrorBar(1:10,nanmean(squeeze(inweight(:,1,:)),2),sem(squeeze(inweight(:,1,:)),2),'lineprops',{'color','r'});
shadedErrorBar(1:10,nanmean(squeeze(inweight(:,2,:)),2),sem(squeeze(inweight(:,2,:)),2),'lineprops',{'color','b'});
shadedErrorBar(1:10,nanmean(squeeze(inweight(:,3,:)),2),sem(squeeze(inweight(:,3,:)),2),'lineprops',{'color','k'});
title('mixing per ordered pred dim (b=multi,r=mixed,k=solo)','fontweight','normal'); 
xlabel('predictive dimension');
ylabel('percent of each type (avg across regions')






















%%
%Camden - timeless
%@input folder is the folder containing the analyzed CCA
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
%grab data
rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_mo\w*.mat'],0,{folder}); 
data_rrr = cellfun(@(x) load(x),fn);

% rec_name = '';
% [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_MSE\w*.mat'],0,{folder}); 
% data_rrr = cellfun(@(x) load(x),fn);

rec_name = '';
[fn,~] = GrabFiles([rec_name '\w*CCA_muaflag1\w*.mat'],0,{folder}); 
data_cca = cellfun(@(x) load(x,'r','paired_areas'),fn);

rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GR\w*.mat'],0,{folder}); 
data_rrr = cellfun(@(x) load(x),fn);

data_rrr(2)=[];
data_cca(2)=[];
%% main
%remove motif 2 (noise)
[area_label,d_opt,d,n] = localDimensionality(data_rrr); 

plot_localDimensionality(area_label,d,n,0); sgtitle('local dimensionality (full)')
plot_localDimensionality(area_label,d_opt,n,1); sgtitle('local dimensionality (optimal)')

%get subspaces and stats
[good_idx,stats] = isSubspace(data_rrr, 0);
% [good_idx,stats] = isSubspace(data_rrr(end), 1);

%Subspace dimensionality
x = MaskSubSpaceStat(cat(2,stats.ss_dim),good_idx,0);
fh = visualizeStrengthMap(x,data_rrr,[0 5]);

%Subspace fraction of full model explained
x = MaskSubSpaceStat(cat(2,stats.prct_pev),good_idx,0);
fh = visualizeStrengthMap(x,data_rrr,[0.7 1]);

%Subspace pctvar
x = cellfun(@nanmean, cat(2,stats.rrr_score));
x = MaskSubSpaceStat(x,good_idx,0);
fh = visualizeStrengthMap(x,data_rrr);

%%
rec_name = LoadDataDirectories(1);
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspacesTemp';
% CreateExampleSubspaceFigures(data_rrr(2),data_cca(2),savedir, rec_name)
CreateExampleSubspaceFigures(data_rrr(end),data_cca(end),savedir, rec_name)

%% the full recordings
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
%grab data
rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GR\w*.mat'],0,{folder}); 
data_rrr = cellfun(@(x) load(x),fn);
data_rrr(2)=[];

%% plot the beta weights for one region
%compare the explained variance
y = arrayfun(@(n) 1-nanmean(data_rrr(n).cvl_rrr{1}), 1:size(data_rrr,2),'UniformOutput',0);
e = arrayfun(@(n) sem(data_rrr(n).cvl_rrr{1}), 1:size(data_rrr,2),'UniformOutput',0);
figure; hold on; 
cellfun(@(a,b) errorbar(1:numel(a),a,b,'o--'), y,e)

%optimal number of dimensions
inweight=[];
dss=[];
for j=1:8
    d = arrayfun(@(n) ModelSelect([nanmean(data_rrr(n).cvl_rrr{j}); std(data_rrr(n).cvl_rrr{j})/sqrt(size(data_rrr(n).cvl_rrr{j},1))], 1:size(data_rrr(n).cvl_rrr{j},2)), 1:size(data_rrr,2),'UniformOutput',1);
    dss{j} = d;
    B = arrayfun(@(n) data_rrr(n).rrr_B{j}(:,1:d(n)), 1:size(data_rrr,2),'UniformOutput',0);

    %plot the first dimension for each one
%     close all; 
    % for dd = 1:min(d)
%     dd=1;
%     figure('units','normalized','position',[0 0 1 1]); hold on; 
%     for i = 1:numel(B)
%         subplot(4,4,i); hold on; 
%         bar(B{i}(:,dd));        
%     end
%     sgtitle(sprintf('dimension %d',dd));
    % end

    grp = data_rrr(1).grouping{j};
    [d,contrib,p] = cellfun(@(x) BetaMixing(x,grp),B,'UniformOutput',0);
    inweight{j} = cat(1,d{:});
end
%get mixing
inweight = cat(3,inweight{:});
a = squeeze((inweight(:,3,:)))'
b = squeeze((inweight(:,2,:)))'
c = squeeze((inweight(:,1,:)))'
dss = cat(1,dss{:});

[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_mo\w*.mat'],0,{folder}); 
data_temp = cellfun(@(x) load(x),fn);
[~,~,d] = localDimensionality(data_temp); 

%plot average dimensionality of a brain region and it's com subspace
d(:,2)=[];
% d=nanmean(d,2)

[rhoa,p] = corr(a(:),d(:))
[rhoa,p] = corr(b(:),d(:))
[rhoa,p] = corr(c(:),d(:))
[rhoa,p] = corr(dss(:),d(:))

%% get the strength mappings
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
%grab data
rec_name = {'Mouse 331 Recording 1','Mouse 331 Recording 2','Mouse 332 Recording 1','Mouse 332 Recording 2','Mouse 334 Recording 1'};


[fn,~] = GrabFiles([rec_name{1} '\w*RRR_muaflag1_GR\w*.mat'],0,{folder}); 
data_rrr = cellfun(@(x) load(x,'contribution','area_label','paired_areas','inactive_idx'),fn(1));
TruePairings = data_rrr(1).paired_areas;
TrueArea = data_rrr(1).area_label;

s_all=[];
for i = 1:numel(rec_name)
    [fn,~] = GrabFiles([rec_name{i} '\w*RRR_muaflag1_GR\w*.mat'],0,{folder}); 
    data_rrr = cellfun(@(x) load(x,'contribution','area_label','paired_areas','inactive_idx'),fn);
    data_rrr(2)=[];
    if size(data_rrr,2)==15
        data_rrr(end)=[];
    end
    xPairings = data_rrr(1).paired_areas;
    xArea = data_rrr(1).area_label;
    
    s = arrayfun(@(n) data_rrr(n).contribution, 1:size(data_rrr,2),'UniformOutput',0);
    s = cellfun(@(x) cat(1,x{:}), s,'UniformOutput',0);
    s = cellfun(@(x) x(:),s,'UniformOutput',0);
    s = cellfun(@(x) MatchToArea(TruePairings,TrueArea,xPairings, xArea,x),s,'UniformOutput',0);
    s = cat(2,s{:});    
    s_all(:,:,i)=s;
%     s = zscore(s,0,3);
%     figure; hold on;
%     for j = 1:size(s,3)
%        subplot(4,4,j);
%        imagesc(s(:,:,j)); 
%        colorbar
%     end
end

%loop through each one
p=NaN(size(s_all,1),1);
for i = 1:size(s_all,1)
   temp = squeeze(s_all(i,:,:)); 
   temp(:,isnan(temp(1,:)))=[]; 
   temp = zscore(temp,0,1);
   grp = repmat([1:size(temp,1)]',1,size(temp,2));
   [p(i),~,stats]=anova1(temp(:),grp(:),'off');
%    imagesc(temp);
%    pause();
end
%% anova testing for differences in strength
p=[];
for j = 1:8
s = arrayfun(@(n) data_rrr(n).cvl_rrr_unique{1}, 1:size(data_rrr,2),'UniformOutput',0);
s = cat(3,s{:});
%compare each strength
for i = 1:7
   temp = squeeze(s(:,i,:));
   [p(i,j),~,stats]=anova1(temp,[],'off');
%    c = multcompare(stats);   
end
end








%

% tss = [];
% for i = 1:size(paired_areas,1)
%     fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,1));
%     x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
%     y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};
% 
%     %normalize to baseline
%     x = normalizeToBaseline(x,[1:2],'mean');
%     y = normalizeToBaseline(y,[1:2],'mean');
% 
%     %use post stimulus
%     x = x(:,3:end,:);
%     y = y(:,3:end,:);
% 
%     %subtract the psth
%     x = x-nanmean(x,3);
%     y = y-nanmean(y,3);
% 
%     %concatentate across trials and pca
%     x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
%     y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
% 
%     N = size(y, 1);
%     tss(i) = sum( sum( ( y - repmat( mean(y), [N 1] ) ).^2 ) );
% end %subspace identification loop

[~,n] = numActiveNeu(data_rrr);
figure; plot(n(:,1),tss,'o')


figure; hold on; 
x = MaskSubSpaceStat(cat(2,stats.prct_pev),good_idx,0);
lm = fitlm(x(:),n(:));
plot(lm)

x = cellfun(@nanmean, cat(2,stats.rrr_score));
x = MaskSubSpaceStat(x,good_idx,0);
lm = fitlm(x(:),n(:));
plot(lm)


























%%
% %get the strength
% [cv_map,x] = cvStrengthMap(data,'r_first');
% for i = 
% 
% 


%
%for each subspace, grab the strongest motif and compare to the 
% 
% %Area list
% area_label = data(1).area_label;
% paired_areas = data(1).paired_areas;
% [cort_ss,cort_idx] = GetCorticalSubspaces(paired_areas, area_label);
% 
% a = [data.a]; 
% idx = cellfun(@(x) isempty(x),a);
% a(idx) = {NaN};
% a = cellfun(@(x) x(:,1),a,'UniformOutput',0);
% %loop through each subspace and compare the angle between coeficients
% for i = 1:size(cort_idx,1)   
%    temp = a(cort_idx(i),:);   
%    %remove empty ones
%    temp(cellfun(@(x) sum(isnan(x))==1,temp,'UniformOutput',1))=[];
%    temp = [temp{:}];
%    cstat = SubspaceConsistency(temp,temp,'angle',1);
%    rng('default')
%    cstat_perm = NaN(1,1000);
%    for j = 1:1000
%       perm_idx = arrayfun(@(n) randperm(size(temp,1),size(temp,1)),1:size(temp,2),'UniformOutput',0);
%       perm_idx = cat(1,perm_idx{:})';
%       cstat_perm(j) = nanmean(SubspaceConsistency(temp(perm_idx),temp,'angle',1));
%    end
% end

%to comparisons. One would be within similarity at different timepoints
%the other would be more othaganol then expected by chance


%things that are important
%proper normalization to baseline
%per timepoint
%remove diagnol
%probably should do better grouping by brain areas

%compare the average r squares or the median (strength)
%compare the number (also strength)
%Diagnols are odd

%other things to plot


% %note camden, this is doing the cca per timepoint
% %get the similarity in angle between all significant timepoints
% %(variability suggests not a valid cc)
% %show a heatmap of when the significant cvs occur
% idx = data(1).t{1};
% for i = 1:numel(data)  
%     %get strength
%     r = data(i).r;
%     temp_idx = cellfun(@(x) isempty(x),r);
%     r(temp_idx) = [];  
%     r = cellfun(@(x) x(1),r,'UniformOutput',1);  
%     r = repmat(r,1,2)';
%     r = r(:);
% 
%     %get number per time
%     temp = data(i).best_idx;
%     temp_idx = cellfun(@(x) isempty(x),temp);
%     temp(temp_idx) = [];  
%     temp = cellfun(@(x) x(1),temp,'UniformOutput',1);
%     
%     temp = arrayfun(@(x) idx(:,x),temp,'UniformOutput',0);
%     
%     temp = cat(1,temp{:});
%     
%     figure; 
%     plot(temp,r,'linestyle','none','marker','x')
%     title(sprintf('%d',i));
% end
% 
% %compare the strength of the best CC for each motif 
% [cv_map,x] = cvStrengthMap(data,'r_first');
% fh = visualizeStrengthMap(cv_map,data);
% 
% cv_map = cvStrengthMap(data,'r_all_sum');
% fh = visualizeStrengthMap(cv_map,data,[]);
% 
% %similarity between upper and lower triangles


%get indices in upper triangle
%get indices in lower triangle
%compare the correlation coefs of the two




















% 
% motif_cvs = cat(2,motif_cvs.motif_cvs);
% cv_map = cvStrengthMap(data,motif_cvs,'r_norm');
% fh = visualizeStrengthMap(cv_map,data,[]);

% %visualize the strength of subspaces per motif 
% cv_map = cvStrengthMap(data,'r');
% % cv_map = cv_map-repmat(nanmean(cv_map,3),1,1,size(cv_map,3));
% % fh = visualizeStrengthMap(cv_map,data,[-0.05 0.05]);
% fh = visualizeStrengthMap(cv_map,data,[0 0.5]);
% 
% cv_map = cvStrengthMap(data,'r_norm');
% fh = visualizeStrengthMap(cv_map,data,[0 0.5]);
% 
% %visualize the significance
% cv_map = cvStrengthMap(data,'pval');
% fh = visualizeStrengthMap(-log(cv_map),data,[]);
% 
% cv_map = cvStrengthMap(data,'num_neu');
% fh = visualizeStrengthMap(cv_map,data,[]);
% 
% 
% cv_map = cvStrengthMap(data,'cv_weight_norm');
% fh = visualizeStrengthMap(cv_map,data,[0 0.5]);

% Question 2: do they follow patterns of activity that one might expect


%compare the peak to the r



%first test is to ask the question: are there significant subspaces and
%does the strength of a subspace match what we would expect given their
%motifs

%take one example pair of regions and test whether it is weaker in motifs
%that do not predominately egage that region. 

%could then be done at a larger scale by clustering the subspace maps
%(correlation) and showing that motifs with similar strength subspaces have
%similar spatial patterns


%the next question to ask is, within a motif, is it a global 


%% 
%plot the average trajectories of each subspace
% 
% 
% %map of significant subspaces
% sub_map = NaN(size(area_label,2));
% for i = 1:size(paired_areas,1)
%     sub_map(paired_areas(i,1),paired_areas(i,2))=r{i}(1);
% end
% close all; 
% figure; 
% imagesc(sub_map,[0 1]); colormap magma; colorbar
% set(gca,'ytick',1:numel(area_label),'YTickLabel',area_label);
% set(gca,'xtick',1:numel(area_label),'XTickLabel',area_label,'XTickLabelRotation',45);
% 
% 
% %plot the variance explained across populations versus the pca angle
% figure; hold on; 
% temp = cellfun(@(x) x(1),r);
% lm = fitlm(cat(1,pca_theta(:,1),pca_theta(:,2)),cat(1,temp,temp));
% plot(lm)
% 
% % [stat,stat_mat,g] = CompareSubspacesWithinAPopulation(cat(1,a,b),cat(1,paired_areas(:,1),paired_areas(:,2)),'angle');
% % temp = cat(1,stat{:});
% % close all; 
% % figure; 
% % imagesc(temp,[0 90]); colormap magma; colorbar
% % set(gca,'ytick',1:numel(g),'YTickLabel',area_label(g));