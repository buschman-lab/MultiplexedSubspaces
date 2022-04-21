function Plot_CompareBetas(data,folder)
%Camden MacDowell - timeless

folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\SubspaceComparison';
beta = cell(1,6);
for i = 1:6
    [fn,~] = GrabFiles(['\w*shuf\w*_rec',num2str(i),'.mat'],0,{folder}); 
    beta{i} = cellfun(@(x) load(x,'rsq','V','B'),fn);
end
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';

data = cell(1,6);
for cur_rec = 1:6
    rec_name = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDm\w*.mat'],0,{folder}); 
    data{cur_rec} = cellfun(@(x) load(x,'cvl_ridge','cvl_rrr','area_label','motif','rrr_V','rrr_B','area_val'),fn);
end

%remove noise motif and null motif
for i = 1:6
   idx = ismember(arrayfun(@(n) data{i}(n).motif, 1:size(data{i},2)),[2,16]);
   data{i}(idx)=[];     
end

%% adjust for the fact no RSP and bfd in rec 3 and 4
for i = 1:size(beta{3},2)
   temp = NaN(8,14,14,10);
   temp([1:3,5,7,8],:,:,:) = beta{3}(i).rsq;
   beta{3}(i).rsq = temp;
end
for i = 1:size(beta{4},2)
   temp = NaN(8,14,14,10);
   temp([1:3,5:8],:,:,:) = beta{4}(i).rsq;
   beta{4}(i).rsq = temp;
end

%% 





















%RR
[~,area_label] = LoadVariable(data,'rrr_dim',[]);

trial_var = LoadVariable(data,'ttt_variability',[]);
% trial_act = LoadVariable(data,'engagement',[]);

for cur_a = 1:numel(area_label)   
    tempvar = squeeze(trial_var(:,:,strcmp(area_label,area_label{cur_a})==1));
%     tempact = squeeze(trial_act(:,:,strcmp(area_label,area_label{cur_a})==1));
    for cur_d = [1,2,3,5]
        figure; hold on;
        t=tiledlayout(4,3); t.TileSpacing = 'compact'; t.Padding = 'compact';    

        %get predictive power across motifs 
        rsq = LoadVariable(data,'ridge_performance',area_label{cur_a});              

        %Split data into two groups
        rng('default')
        clustID = kmeans(cat(1,nanmean(rsq),nanmean(tempvar))',2);
        
        %within each group, sort by increasing weighted rsq        
        rsq = rsq.*tempvar;
        [~,idx_grpone] = sort(nanmean(rsq(:,clustID==1)),'ascend');
        [~,idx_grptwo] = sort(nanmean(rsq(:,clustID==2)),'ascend');
        grpone = find(clustID==1);
        grptwo = find(clustID==2);
        strongidx=[];
        if nanmean(rsq(grpone(:)))>nanmean(rsq(grptwo(:)))  %one > two        
            idx = [grptwo(idx_grptwo)',grpone(idx_grpone)'];
            strongidx(1,:) = idx(numel(grptwo)+1:end);
            strongidx(2,:) = numel(grptwo)+1:numel(idx);
        else
            idx = [grpone(idx_grpone)',grptwo(idx_grptwo)'];
            strongidx(1,:) = idx(numel(grpone)+1:end);
            strongidx(2,:) = numel(grpone)+1:numel(idx);
        end
        

        nexttile([1,3])
        plot(nanmean(rsq(:,idx)))
        fp.FormatAxes(gca);
        title(sprintf('Contribution of dimension %d',cur_d),'fontweight','normal')
        set(gca,'xtick',1:14,'xticklabel',idx(1:14))
        xlabel('motif'); ylabel({'weighted','performance'});        
        
        %get the rsq for brain area        
        xbeta = arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
        xbeta = nanmean(cat(3,xbeta{:}),3);
        xbeta = xbeta(idx,idx);

        %the top value from each permuted distribution
        permdist = 1:1000;
        for cur_perm = 1:1000
            perm = arrayfun(@(n) squeeze(beta{n}(cur_perm).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
            permdist(cur_perm) = max(nanmean(cat(3,perm{:}),3),[],'all');
        end
        alphamap = double(xbeta>prctile(permdist,95));
        alphamap(alphamap==0)=0.15;
        
        nexttile([3,3])
        h = imagesc(xbeta); c=colorbar; colormap(magma)
        ylabel(c,'percent explained variance');
        h.AlphaData = alphamap;
        set(gca,'xtick',1:numel(idx),'xticklabel',idx,'ytick',1:numel(idx),'YTickLabel',idx,'ydir','normal')
        xlabel('motif'); ylabel('motif')
        set(gcf,'units','centimeters','position',[8 8 12 12])
        title(t,sprintf('Area %s',area_label{cur_a}));        
        %add a box around 'strong' and weak
        rectangle('position',[strongidx(2,1)-0.5,strongidx(2,1)-0.5,size(strongidx,2)+1,size(strongidx,2)+1],'LineStyle',':','EdgeColor',fp.c_lr,'LineWidth',1.5)
        rectangle('position',[0.5 0.5,strongidx(2,1)-1,strongidx(2,1)-1],'LineStyle',':','EdgeColor',fp.c_glm,'LineWidth',1.5)        
    end
end


handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng'},'example figures',savedir,0); %close all 


%% Per recording and then averaged (this is the correct way to do it, fyi)
%NEED TO MAKE SURE we are averging by the motif that is not being change
%(i.e. column in rsq)
sigval = cell(6,1);
strongPEV = cell(6,numel(area_label));
weakPEV = cell(6,numel(area_label));
rrr_strong = cell(6,numel(area_label));
rrr_weak = cell(6,numel(area_label));
strongPEV_col = cell(6,numel(area_label));
for cur_rec = 1:6
for cur_a = 1:numel(area_label)   
    tempvar = squeeze(trial_var(:,:,strcmp(area_label,area_label{cur_a})==1));
    %get predictive power across motifs 
    rsq = LoadVariable(data,'ridge_performance',area_label{cur_a}); 
    rrr_dim = LoadVariable(data,'rel_performance',area_label{cur_a}); 

    %Split data into two groups (consistent across recordings)
    rng('default')
    clustID = kmeans(cat(1,nanmean(rsq),nanmean(tempvar))',2);

    %within each group, sort by increasing weighted rsq        
    rsq = rsq.*tempvar;
    [~,idx_grpone] = sort(nanmean(rsq(:,clustID==1)),'ascend');
    [~,idx_grptwo] = sort(nanmean(rsq(:,clustID==2)),'ascend');
    grpone = find(clustID==1);
    grptwo = find(clustID==2);
    strongidx=[];
    if nanmean(rsq(grpone(:)))>nanmean(rsq(grptwo(:)))  %one > two        
        idx = [grptwo(idx_grptwo)',grpone(idx_grpone)'];
        strongidx(1,:) = idx(numel(grptwo)+1:end);
        strongidx(2,:) = numel(grptwo)+1:numel(idx);
    else
        idx = [grpone(idx_grpone)',grptwo(idx_grptwo)'];
        strongidx(1,:) = idx(numel(grpone)+1:end);
        strongidx(2,:) = numel(grpone)+1:numel(idx);
    end

    temp = squeeze(rrr_dim(cur_rec,:,:));
    rrr_strong{cur_rec,cur_a} = temp(strongidx(1,:),:);
    rrr_weak{cur_rec,cur_a} = temp(strongidx(1,:),:);
    
    for cur_d = 1:10        
        %get the rsq for brain area        
        xbeta = squeeze(beta{cur_rec}(1).rsq(cur_a,:,:,cur_d));
        xbeta = xbeta(idx,idx);

        %the top value from each permuted distribution
        permdist = arrayfun(@(n) max(squeeze(beta{cur_rec}(n).rsq(cur_a,:,:,cur_d)),[],'all'), 1:1000,'UniformOutput',1);
       
%         temp = xbeta(strongidx(2,:),strongidx(2,:));
        temp = xbeta;                
        strongPEV{cur_rec,cur_a}(cur_d,:) = temp(~isnan(temp));
        strongPEV_col{cur_rec,cur_a}(cur_d,:) = nanmean(temp);

        temp = xbeta(1:strongidx(2,1)-1,1:strongidx(2,1)-1);
        weakPEV{cur_rec,cur_a}(cur_d,:) = temp(~isnan(temp));
        
        %strong

        %95% significance value
        sigval{cur_rec}(cur_a,cur_d) = prctile(permdist,95);
    end
    
end
end %rec

%% number of dimensions before dropping off the map
close all;
pevall = NaN(8,1000);
for cur_a = 1:8
    tempall = cell(1,6);
    for cur_rec = 1:6
        temp = strongPEV_col{cur_rec,cur_a}-repmat(sigval{cur_rec}(cur_a,:)',1,size(strongPEV_col{cur_rec,cur_a},2));
        temp = arrayfun(@(n) find(temp(:,n)>0,1,'last'),1:size(temp,2),'UniformOutput',0);
        a = cellfun(@(x) isempty(x),temp);
        temp(a)={0};
        temp = [temp{:}];        
        b = [zeros(numel(temp),1),rrr_strong{cur_rec,cur_a}];
        tempall{cur_rec} = arrayfun(@(n) b(n,temp(n)+1),1:numel(temp),'UniformOutput',1);
    end
    tempall = nanmean(cat(1,tempall{:}),2);
    pevall(cur_a,1:numel(tempall)) = tempall;
end

pevall(:,sum(isnan(pevall),1)==8)=[];

figure; hold on; 
col = getColorPalet(8);
[~,x] = sort(nanmean(pevall,2),'descend');
b=cell(1,8);
for i = 1:8
    temp = pevall(x(i),:);
    temp(isnan(temp))=[];
    b{i}=bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
[pval,tbl,stats] = anova1(pevall(x,:)',[],'off');
hh = line(nan, nan, 'Color', 'none');
legend(hh, {sprintf('fstat=%0.3f \n p=%0.2f',tbl{2,5},pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best')     
ylabel({'PEV of subspace dimensions','that are significantly shared'});
set(gca,'xtick',1:8,'xticklabel',area_label(x))
title({'Brain regions differ in how much subspace is dominated by shared subspaces'},'fontweight','normal')



%so we should scale performance by that value
figure; hold on; 
Plot_CompareValueBetweenMotifs(pevall',x,fp,'right');
ylabel({'Weighted','performance'});
title(sprintf('Predicting %s trial-to-trial \n variability from other areas (scaled)',area_name),'fontweight','normal')




handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng'},'Topology',savedir,0); close all    




%% average across recordings
sigval = NaN(numel(area_label),10);
strongPEV = cell(1,numel(area_label));
weakPEV = cell(1,numel(area_label));
rrr_strong = cell(1,numel(area_label));
rrr_weak = cell(1,numel(area_label));
strongPEV_col = cell(1,numel(area_label));
for cur_a = 1:numel(area_label)   
    tempvar = squeeze(trial_var(:,:,strcmp(area_label,area_label{cur_a})==1));
    %get predictive power across motifs 
    rsq = LoadVariable(data,'ridge_performance',area_label{cur_a}); 
    rrr_dim = LoadVariable(data,'rel_performance',area_label{cur_a}); 

    %Split data into two groups
    rng('default')
    clustID = kmeans(cat(1,nanmean(rsq),nanmean(tempvar))',2);

    %within each group, sort by increasing weighted rsq        
    rsq = rsq.*tempvar;
    [~,idx_grpone] = sort(nanmean(rsq(:,clustID==1)),'ascend');
    [~,idx_grptwo] = sort(nanmean(rsq(:,clustID==2)),'ascend');
    grpone = find(clustID==1);
    grptwo = find(clustID==2);
    strongidx=[];
    if nanmean(rsq(grpone(:)))>nanmean(rsq(grptwo(:)))  %one > two        
        idx = [grptwo(idx_grptwo)',grpone(idx_grpone)'];
        strongidx(1,:) = idx(numel(grptwo)+1:end);
        strongidx(2,:) = numel(grptwo)+1:numel(idx);
    else
        idx = [grpone(idx_grpone)',grptwo(idx_grptwo)'];
        strongidx(1,:) = idx(numel(grpone)+1:end);
        strongidx(2,:) = numel(grpone)+1:numel(idx);
    end

    temp = squeeze(nanmean(rrr_dim,1));
    rrr_strong{cur_a} = temp(strongidx(1,:),:);
    rrr_weak{cur_a} = temp(strongidx(1,:),:);
    
    for cur_d = 1:10
        
        %get the rsq for brain area        
        xbeta = arrayfun(@(n) squeeze(beta{n}(1).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
        xbeta = nanmean(cat(3,xbeta{:}),3);
        xbeta = xbeta(idx,idx);

        %the top value from each permuted distribution
        permdist = 1:1000;
        for cur_perm = 1:1000
            perm = arrayfun(@(n) squeeze(beta{n}(cur_perm).rsq(cur_a,:,:,cur_d)), 1:size(beta,2),'UniformOutput',0);
            permdist(cur_perm) = max(nanmean(cat(3,perm{:}),3),[],'all');
        end
        
        %
        temp = xbeta(strongidx(2,:),strongidx(2,:));
        strongPEV{cur_a}(cur_d,:) = temp(~isnan(temp));
        strongPEV_col{cur_a}(cur_d,:) = nanmean(temp);

        temp = xbeta(1:strongidx(2,1)-1,1:strongidx(2,1)-1);
        weakPEV{cur_a}(cur_d,:) = temp(~isnan(temp));
        
        %strong

        %95% significance value
        sigval(cur_a,cur_d) = prctile(permdist,95);
    end
    
end
%% For each area, plot over time
% weakPEV(weakPEV<0)=NaN;
% strongPEV(strongPEV<0)=NaN;
% sigval(sigval<0)=NaN;

%can also plot the fraction above significance

close all;
for cur_a = 1:numel(area_label)
    figure; hold on; 
    col = getColorPalet(8);
    col = arrayfun(@(n) col(n,:),1:8, 'UniformOutput',0);
    temp = weakPEV{cur_a};
    temp(temp<0)=NaN;
    col = repmat(col,size(temp,2),1)';
    vp = CompareViolins(temp,fp,'col',col(1,:),'connectline',[0.25 0.25 0.25 0.50],'plotspread',0,'divfactor',3);
    plot(1:size(sigval,2),sigval(cur_a,:),'linewidth',2,'color',[0.8 0.1 0.1],'linestyle',':')
    ylabel('Generalization Performance'); xlabel('Subspace Dimension')    
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,sprintf('Weak Area %s',area_label{cur_a}))
    fp.FigureSizing(gcf,[2 2 8 4],[8 8 14 10])
    ylim([0 1])

    figure; hold on; 
    col = getColorPalet(8);
    col = arrayfun(@(n) col(n,:),1:8, 'UniformOutput',0);
    temp = strongPEV{cur_a};
    temp(temp<0)=NaN;
    col = repmat(col,size(temp,2),1)';
    vp = CompareViolins(temp,fp,'col',col(1,:),'connectline',[0.25 0.25 0.25 0.50],'plotspread',0,'divfactor',3);
    plot(1:size(sigval,2),sigval(cur_a,:),'linewidth',2,'color',[0.8 0.1 0.1],'linestyle',':')
    ylabel('Generalization Performance'); xlabel('Subspace Dimension')    
    fp.FormatAxes(gca); grid on
    fp.SetTitle(gca,sprintf('Strong Area %s',area_label{cur_a}))
    fp.FigureSizing(gcf,[2 2 8 4],[8 8 14 10])
    ylim([0 1])
end



handles = findall(groot,'Type','figure');
saveCurFigs(handles,{'-dpng'},'across dimensions',savedir,0); close all    



%% number of dimensions before dropping off the map
close all;
pevall = NaN(8,1000);
for cur_a = 1:8   
    temp = strongPEV_col{cur_a}-repmat(sigval(cur_a,:)',1,size(strongPEV_col{cur_a},2));
    temp = arrayfun(@(n) find(temp(:,n)>0,1,'last'),1:size(temp,2),'UniformOutput',0);
    a = cellfun(@(x) isempty(x),temp);
    temp(a)={0};
    temp = [temp{:}];        
    b = [zeros(numel(temp),1),rrr_strong{cur_a}];
    tempall = arrayfun(@(n) b(n,temp(n)+1),1:numel(temp),'UniformOutput',1);
    pevall(cur_a,1:numel(tempall)) = tempall;
end

figure; hold on; 
col = getColorPalet(8);
[~,x] = sort(nanmean(pevall,2),'descend');
b=cell(1,8);
for i = 1:8
    temp = pevall(x(i),:);
    temp(isnan(temp))=[];
    b{i}=bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
[pval,tbl,stats] = anova1(pevall(x,:)',[],'off');
hh = line(nan, nan, 'Color', 'none');
legend(hh, {sprintf('fstat=%0.3f \n p=%0.2f',tbl{2,5},pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best')     
ylabel({'PEV of subspace dimensions','that are significantly shared'});
set(gca,'xtick',1:8,'xticklabel',area_label(x))
title({'Brain regions differ in how much subspace is dominated by shared subspaces'},'fontweight','normal')

%% Why is it the case that some regions have more mixing... one could be the anatomical diversity to other regions (actually this is counter what we see)

ClustInfo = NaN(8,6,3);
for cur_a = 1:8
    fprintf('working on area %d of %d',cur_a,8);
    rrrB = arrayfun(@(n) LoadVariable(data,'rrr_beta_noflip',area_label{cur_a},n),1:15,'UniformOutput',0);
    rrrB = cat(2,rrrB{:});
    rrrV = arrayfun(@(n) LoadVariable(data,'rrr_V_noflip',area_label{cur_a},n),1:15,'UniformOutput',0);
    rrrV = cat(2,rrrV{:});
    rrrF=[];
    for cur_rec = 1:6
        for i = 1:size(rrrB,2)
            b = squeeze(rrrB(cur_rec,i,:));
            v = squeeze(rrrV(cur_rec,i,:));
            temp = b*v';
            rrrF(cur_rec,i,:) = temp(:);
        end
    end

    rrrB = arrayfun(@(n) LoadVariable(data,'rrr_beta',area_label{cur_a},n),1:15,'UniformOutput',0);
    rrrB = cat(2,rrrB{:});
    rrrV = arrayfun(@(n) LoadVariable(data,'rrr_V',area_label{cur_a},n),1:15,'UniformOutput',0);
    rrrV = cat(2,rrrV{:});

    %we see that the correlation between strength and V similarity is stronger
    %than the correlation between strength and B similarity. Suggesting that
    %one reason for this could be that more mixing is because we have
    %overlapping local region representations. 
    for cur_rec = 1:6       
        rng('default')
        a = squeeze(rrrB(cur_rec,:,:));
        if sum(isnan(a(:)))==numel(a)
        else
            a(:,sum(isnan(a),1)>0)=[];
            y = run_umap(a,'n_neighbors',10,'verbose','none','method','java','cluster_method_2D','DBSCAN','metric','spearman','min_dist',0.2,'minpts',3,'epsilon',0.15,'randomize','false');
            idx = dbscan(y,0.3,3);
    %         figure; gscatter(y(:,1),y(:,2),idx)
            ClustInfo(cur_a,cur_rec,1) = max(idx);

            rng('default')
            a = squeeze(rrrV(cur_rec,:,:));
            a(:,sum(isnan(a),1)>0)=[];
            y = run_umap(a,'n_neighbors',10,'verbose','none','method','java','cluster_method_2D','DBSCAN','metric','spearman','min_dist',0.2,'minpts',3,'epsilon',0.15,'randomize','false');
            idx = dbscan(y,0.3,3);
    %         figure; gscatter(y(:,1),y(:,2),idx)
            ClustInfo(cur_a,cur_rec,2) = max(idx);

            rng('default')
            a = squeeze(rrrF(cur_rec,:,:));
            a(:,sum(isnan(a),1)>0)=[];
            y = run_umap(a,'n_neighbors',10,'verbose','none','method','java','cluster_method_2D','DBSCAN','metric','spearman','min_dist',0.2,'minpts',3,'epsilon',0.15,'randomize','false');
            idx = dbscan(y,0.3,3);
    %         figure; gscatter(y(:,1),y(:,2),idx)
            ClustInfo(cur_a,cur_rec,3) = max(idx);
        end
            
    end
end

%%
figure; hold on; 
for i = 1:8
    temp = squeeze(ClustInfo(i,:,1));
    bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
set(gca,'xtick',1:8,'XTickLabel',area_label);

figure; hold on; 
for i = 1:8
    temp = squeeze(ClustInfo(i,:,2));
    bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
set(gca,'xtick',1:8,'XTickLabel',area_label);

figure; hold on; 
for i = 1:8
    temp = squeeze(ClustInfo(i,:,1)./ClustInfo(i,:,2));
    bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
set(gca,'xtick',1:8,'XTickLabel',area_label);

temp = squeeze(ClustInfo(:,:,1)./ClustInfo(:,:,2))';

[~,p] = ttest(temp(:,1),temp(:,5)),

%%

rsq = cellfun(@(x) x(1).rsq,beta,'UniformOutput',0);
rsq = nanmean(cat(5,rsq{:}),5);
V = cellfun(@(x) x(1).V,beta,'UniformOutput',0);
V = nanmean(cat(5,V{:}),5);
B = cellfun(@(x) x(1).B,beta,'UniformOutput',0);
B = nanmean(cat(5,B{:}),5);

a = reshape(squeeze(rsq(:,:,:,cur_d)),8,14*14);
b = reshape(squeeze(B(:,:,:,cur_d)),8,14*14);
v = reshape(squeeze(V(:,:,:,cur_d)),8,14*14);

%one hypothesis is areas with stronger subspace have overlapping networks
%with the rest of the brain between motifs

%compare the betas across all areas versus comparing
cur_d = 2;
a = reshape(squeeze(rsq(:,:,:,cur_d)),8,14*14);
b = reshape(squeeze(B(:,:,:,cur_d)),8,14*14);
v = reshape(squeeze(V(:,:,:,cur_d)),8,14*14);

figure; hold on; 
col = getColorPalet(8);
x = v;
b=cell(1,8);
for i = 1:8
    temp = x(i,:);
    temp(isnan(temp))=[];
    b{i}=bar(i,nanmean(temp),'facecolor',col(i,:),'facealpha',0.5,'EdgeAlpha',0);
    plot(i*ones(numel(temp),1),temp,'color',col(i,:),'marker','o','markersize',fp.markersizesmall,'linestyle','none');
end
[pval,tbl,stats] = anova1(x',[],'off');
hh = line(nan, nan, 'Color', 'none');
legend(hh, {sprintf('fstat=%0.3f \n p=%0.2f',tbl{2,5},pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best')     
ylabel({'PEV of subspace dimensions','that are significantly shared'});
set(gca,'xtick',1:8,'xticklabel',area_label)
title({'Brain regions differ in how much subspace is dominated by shared subspaces'},'fontweight','normal')




%the other hypothesis is that 


%an alternative would be both

cur_a = 5;
m = [3,6]; 

cur_d = 2;
%what is their PEV across dimensions
a = rsq(cur_a,3,6,cur_d);
b = B(cur_a,3,6,cur_d);
v = V(cur_a,3,6,cur_d);

figure; hold on; 
cur_a = 5;
a = squeeze(rsq(cur_a,:,:,cur_d));
b = squeeze(B(cur_a,:,:,cur_d));
v = squeeze(V(cur_a,:,:,cur_d));
plot(fitlm(a(:),v(:)));

figure; hold on; 
plot(fitlm(a(:),b(:)));

%across all pairings, compare correlation in Beta and correlation in V to
%sharedness
a = rsq(8,5,7,1);
b = B(8,5,7,1);
v = V(8,5,7,1);
















end %function end
