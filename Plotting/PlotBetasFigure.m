function PlotBetasFigure(data,savedir)
%Camden MacDowell - timeless
fp = fig_params_cortdynamics;
if nargin <2; savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BetaWeights'; end 
if ~exist(savedir,'dir'); mkdir(savedir); end
%%
fp = fig_params_cortdynamics;
close all;
targ_area='VIS';
cur_d = 1;
figure; hold on; 
cur_rec =5;
cur_m = 7;
[b,area,area_label] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area);
bar(b,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'BarWidth',1,'FaceAlpha',1)
plot(1:numel(b),b,'linewidth',1.5,'color','r')
plot([-20,numel(b)+20],[0,0],'linewidth',1.5,'color','k')
fp.FormatAxes(gca); box on; grid on; 
set(gca,'ydir','reverse','YAxisLocation','right')
xlabel('neurons');
ylabel('subspace betas');
title(sprintf('VIS rec %d m %d',cur_rec,cur_m),'fontweight','normal');
xlim([1,numel(b)])
% set(gca,'ylim',[-0.05 0.15],'ytick',[-0.05:0.05:0.15])
camroll(-90)
fp.FigureSizing(gcf,[3 2 3 4],[2 10 15 15])

%plot the identify on top
id = unique(area);
for i = 1:7
    figure; hold on; 
    [y,edges] = histcounts(find(area==id(i)),0:20:numel(area)+20);   
    %interpolate
    y = interp1(edges(2:end),y,1:numel(area));
    y = smoothdata(y,'gaussian',50);
    %append the last data point for plotting
    y(isnan(y))=0;
    y = [0,y,0];
    y = 100*(y/nansum(y(:))); %normalize to % 
    patch(1:numel(y),y,fp.c_area(id(i),:),'edgecolor','none','facealpha',0.5)
    ylim([0 0.25]);
    set(gca,'clipping','off','ytick',[0 0.25])
    fp.FormatAxes(gca); box off
    xlabel('neurons');
    ylabel('% neurons');
    title(sprintf('%s dim %d r%d m%d',area_label{i},cur_d,cur_rec,cur_m),'fontweight','normal');    
%     xlim([-20,numel(b)+20])
    xlim([1,numel(b)]) 
    camroll(-90)
    fp.FigureSizing(gcf,[3 2 0.75 4],[2 10 15 15])
    set(gca,'YAxisLocation','right')
end
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('NeuronContributions_Prediction%s_dim%d',targ_area,cur_d),savedir,0); close all

%% Get the fraction of top beta weights contributed by each area across motifs and recordings
preload data
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
ent = [];
for cur_d = 1:3
    for cur_area = 1:numel(area_all)
        targ_area= area_all{cur_area};
        B = cell(6,14);
        AArea = cell(6,14);
        AreaLabel = cell(6,14);

        for cur_rec = 1:6
            for cur_m = 1:14        
                [B{cur_rec,cur_m},AArea{cur_rec,cur_m},AreaLabel{cur_rec,cur_m}] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area);
            end
        end
        
        %sweep fractions
        idx = 0:0.1:1;
        x = NaN(7,numel(idx));
        ci = NaN(7,2,numel(idx));
        y = NaN(1,numel(idx));
        xboot = NaN(7,1000,11);
        for i = 1:numel(idx)
            [~,~,stats,y(i)] = SubspaceUniformity(B,AArea,AreaLabel,data,targ_area,idx(i),1);
            x(:,i) = stats.mean;
            ci(:,:,i) = stats.ci';
            xboot(:,:,i) = stats.xboot';            
            %get the average entropy
        end
        ent(cur_d,cur_area,1:numel(idx)) = y;

        %plot the chart
        figure; hold on; 
        col = fp.c_area; col = col(strcmp(area_all,area_all{cur_area})==0,:);
        for i = 1:7    
            tempci = abs(flipud(squeeze(ci(i,:,:))-x(i,:)));
            shadedErrorBar(idx,x(i,:),tempci,'lineprops',{'color',[col(i,:) 0.75],'linewidth',2},'patchSaturation',0.15);
        end
        ylabel('Percentage of Population');
        xlabel('Fraction of Beta weights');
        fp.FormatAxes(gca); box on; grid on
        fp.FigureSizing(gcf,[3 3 2.9 4],[10 10 12 10])
        title(sprintf('%s dimesion %d',targ_area,cur_d));    
        
        %plot the AUC in a small panel
        %mean auc
        xx = x/100;
        auc = NaN(1,size(x,1));
        aucci = NaN(2,size(x,1));
        aucboot = NaN(1000,size(x,1));
        
        for i = 1:size(x,1)
            auc(i) = trapz(xx(i,:)/numel(xx(i,:)));
            temp = squeeze(ci(i,1,:))/100;
            aucci(1,i) = trapz(temp/numel(temp));
            temp = squeeze(ci(i,2,:))/100;
            aucci(2,i) = trapz(temp/numel(temp));
            temp = squeeze(xboot(i,:,:));
            temp = arrayfun(@(n) trapz((temp(n,:)/100)/numel(temp(n,:))),1:size(temp,1),'UniformOutput',1);
            aucboot(:,i) = temp;
        end
        
        %mini plot        
        a = area_all(strcmp(area_all,area_all{cur_area})==0);
        [auc,idx] = sort(auc,'descend');
        a = a(idx);
        col = fp.c_area; col = col(strcmp(area_all,area_all{cur_area})==0,:);
        col = col(idx,:);
        aucci = aucci(:,idx);
        figure; hold on; 
        for i = 1:numel(a)
            errorbar(i,auc(i),auc(i)-aucci(1,i),aucci(2,i)-auc(i),'linestyle','none','marker','o',...
                'MarkerEdgeColor','none','MarkerFaceColor',col(i,:),'markersize',2,'linewidth',1.5,'color',col(i,:),'CapSize',4)            
        end
        xlim([0 8])
        plot([0,8],[0.5 0.5],'linestyle','--','color',[0.25 0.25 0.25],'linewidth',1);
        ylim([0.35,0.6])
        fp.FormatAxes(gca); box on; grid on
        fp.FigureSizing(gcf,[3 3 1.5 1.5],[10 10 12 10]) 
        
        aucboot = aucboot(:,idx);
        %do significance testing: 
        pval = NaN(7,7);        
        figure; hold on;        
        for i = 1:size(pval,1)
            for j= 1:size(pval,2)
               temp = aucboot(:,i)-aucboot(:,j);
               pval(i,j) = sum(temp>0)/numel(temp(~isnan(temp)));
               plot(i,j,'marker','.');
               text(i,j,sprintf('%0.3f',pval(i,j)),'FontWeight','normal','fontsize',10);
            end
        end
        
        
        
    end
    
    e = squeeze(ent(cur_d,:,:));
    figure; hold on; 
    col = fp.c_area; 
    p = repmat(idx,7,1);
    nullent = arrayfun(@(n) -sum(p(:,n).*log2(p(:,n))),1:size(p,2)); 
    plot(idx,ones(1,numel(idx)),'linestyle',':','color',[0.75 0.75 0.75],'linewidth',2)
    for i = 1:8
        plot(idx,e(i,:)./nullent,'color',col(i,:),'linewidth',2)
    end
    xlim([idx(2),idx(end-1)]);
    ylabel('Relative Entropy');
    xlabel('Fraction of Beta weights');
    fp.FormatAxes(gca); box on; grid on
    fp.FigureSizing(gcf,[3 3 2 3.15],[10 10 12 10])
    title(sprintf('Uniformity of incoming \ninteractions Dimesion %d',cur_d));     

end

%
set(gca,'ylim',[0 50],'ytick',[0,25,50])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('AreaContributions_Prediction%s_dim%d',targ_area,cur_d),savedir,0); close all
save([savedir,filesep,sprintf('AreaContributions_Prediction%s_dim%d',targ_area,cur_d)],'stats','n');

%% Plot the correlation between 
% do two examples. One, Everything predicting MOs or VIS another everything predicting PRE or RSP and compare the
% betas within each region when predicting the other. Plot as a correlation
% chart
close all;
fp = fig_params_cortdynamics;
%get beta in all areas predicting 
targ_areas = {'PRE','VIS'};
cur_d = 1;

cur_rec = 2; 
cur_m = 1; %def: 5
% Plot and example: HIPP or RSP to VIS and SS boostrap confidence intervals and significance
ExampleSubspaceGen(data,cur_rec,cur_m,cur_d,'HIPP',targ_areas)
% ExampleSubspaceGen(data,5,1,cur_d,'THAL',targ_areas)
ExampleSubspaceGen(data,5,7,3,'MOs',{'VIS','SS'})
targ_areas = {'VIS','SS'};
ExampleSubspaceGen(data,cur_rec,cur_m,[2,4],'MOs',targ_areas)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('ExampleProjections_%sand%s_dim%d',targ_areas{1},targ_areas{2},cur_d),savedir,0); close all

%% plot examples of ones that are not particularly well shared
% targ_areas = {'PRE','VIS'};
% ExampleSubspaceGen(data,cur_rec,cur_m,[3,2],'HIPP',targ_areas)
targ_areas = {'VIS','SS'};
ExampleSubspaceGen(data,cur_rec,cur_m,[3,2],'MOs',targ_areas)
fp.FigureSizing(gcf,[3 2 3 4],[2 10 15 15])
% for i = 1:4
%     for j= 1:4
%         ExampleSubspaceGen(data,cur_rec,cur_m,[i,j],'MOs',targ_areas)
%     end
% end
%%
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('ExampleNonShared%sand%s',targ_areas{1},targ_areas{2}),savedir,0); close all

%% run across all recs, motifs, and areas
close all; 
targ_areas = {'VIS','SS'};
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
for cur_d = 1:2
    [rho,rho_perm,~,~,n,nn] = SubspaceSimilarity(data,cur_d,targ_areas);
    [rho,stats] = pairedBootstrap(abs(rho),@nanmean);
    [rho_perm,permstats] = pairedBootstrap(abs(rho_perm),@nanmean);
    % plot both the absolute value and tru value... 
    ContributionPlot(rho,intersect(n,nn),data,[0.2,0.65],rho_perm)
    ylabel('Beta rho_z');
    temp = get(gca,'ylim');
    set(gca,'ylim',[0,temp(2)]);
    title(sprintf('Predicting %s %s DIM %d',targ_areas{1},targ_areas{2},cur_d),'FontWeight','normal');
    saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('CombinedSimilarityPlot%sand%s_dim%d',targ_areas{1},targ_areas{2},cur_d),savedir,0); close all
    save([savedir,filesep,sprintf('CombinedSimilarityPlot%sand%s_dim%d.mat',targ_areas{1},targ_areas{2},cur_d)],'stats','permstats');
end
% do for a different example area
close all
targ_areas = {'PRE','VIS'};
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
for cur_d = 1:2
    [rho,rho_perm,~,~,n,nn] = SubspaceSimilarity(data,cur_d,targ_areas);
    [rho,stats] = pairedBootstrap(abs(rho),@nanmean);
    [rho_perm,permstats] = pairedBootstrap(abs(rho_perm),@nanmean);
    % plot both the absolute value and tru value... 
    ContributionPlot(rho,intersect(n,nn),data,[0.2,0.65],rho_perm)
    ylabel('Beta rho_z');
    set(gca,'ylim',[0,temp(2)]);
    title(sprintf('Predicting %s %s DIM %d',targ_areas{1},targ_areas{2},cur_d),'FontWeight','normal');
    saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('CombinedSimilarityPlot%sand%s_dim%d',targ_areas{1},targ_areas{2},cur_d),savedir,0); close all
    save([savedir,filesep,sprintf('CombinedSimilarityPlot%sand%s_dim%d.mat',targ_areas{1},targ_areas{2},cur_d)],'stats','permstats');
end

%% Plot similarity across all areas per dimension
for cur_d = 1:10
    r = SubspaceSimilarityFull(data,cur_d);
    %equate sizes
    m = max(cellfun(@numel,r));
    rho = cellfun(@(x) cat(1,x,NaN(m-numel(x),1)),r,'UniformOutput',0);
    rho = abs(cat(2,rho{:}));
    rho_all{cur_d} = rho;
    [rho,~] = pairedBootstrap(fisherZ(rho),@nanmean);
    [~, area_label] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
    ContributionPlot(rho,area_label,data,[.2,0.65])
    ylabel('Beta Rho_z')
    title(sprintf('Dimension %d',cur_d),'fontweight','normal');
    
    %generate a graph
    
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'BetaSimilarityAcrossDimensions',savedir,0); close all

%% Plot the average
rho = cat(3,rho_all{:});
rho = rho(:,:,1:5);
rho = fisherZ(rho);
% rho = permute(rho,[2,1,3]);
% rho = reshape(rho,size(rho,1),size(rho,2)*size(rho,3))';
rho = nanmean(rho,3);
[rho,~] = pairedBootstrap(rho,@nanmean);
ContributionPlot(rho,area_label,data,[.2,0.65])
set(gca,'ylim',[0.27,0.43])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'CombinedBetaSimilarityAcrossDimensions',savedir,0); close all


%%
fp = fig_params_cortdynamics;
[~, area_name] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
figure; hold on;       
col = fp.c_area(ismember(area_name,area_label),:);
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
%sort
[~,ind] = sort(nanmean(nneu),'descend');
col = col(ind);
CompareViolins(nneu(:,ind)',fp,'plotspread',0,'divfactor',spreadfactor(1),'plotaverage',1,'col',col,'distWidth',spreadfactor(2));
% if ~isempty(xperm)
%    CompareViolins(xperm(:,ind)',fp,'plotspread',0,'divfactor',spreadfactor(1),'plotaverage',0,'col',repmat({[0 0 0]},1,numel(ind)),'distWidth',spreadfactor(2)); 
% end
set(gca,'XTickLabel',area_label(ind),'XTickLabelRotation',45)
fp.FormatAxes(gca);  box on; grid on; 
fp.FigureSizing(gcf,[3 3 6 3],[10 10 14 10])    
ylabel('% of population');




%% Plot similar network engagement of two areas (based off shared representation in third area)
for cur_d = 1:5
   rmat = SubspaceSimilarityMatrix(data,cur_d);
   %average per rec and across motifs
   r = squeeze(nanmedian(nanmedian(abs(rmat),1),2));
   %average across areas 
   rfull = nanmedian(r,3);
   %symmetric
   rfull = tril(rfull.') + triu(rfull);
   rfull(isnan(rfull))=0;
   close all; figure; hold on; 
   circularGraph(rfull.^2,'colormap',repmat([0.5 0.5 0.5],8,1),'label',data{1}(1).area_label,'normVal',0.15);
   title(sprintf('Dim %d Subspace Networks (rsq range: %0.2f-%0.2f)',cur_d,min(rfull(rfull>0)).^2,max(rfull(:)).^2),'FontWeight','normal')
   set(findall(gca, 'type', 'text'), 'visible', 'on')
   saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('SubspaceNetworks_Combo_dim%d',cur_d),savedir,0); close all
   % do for each region
   for i = 1:8
       %average across areas 
       rfull = r(:,:,i);
       %symmetric
       rfull = tril(rfull.') + triu(rfull);
       rfull(isnan(rfull))=0;
       figure; hold on; 
       circularGraph(rfull.^2,'colormap',repmat(fp.c_area(i,:),8,1),'label',data{1}(1).area_label,'normVal',0.25);
       fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
       saveCurFigs(get(groot, 'Children'),{'-dpng'},sprintf('SubspaceNetworks_Indi_%sdim%d',data{1}(1).area_label{i},cur_d),savedir,0); 
       title(sprintf('Dim %d Subspace Networks (rsq range: %0.2f-%0.2f)',cur_d,min(rfull(rfull>0)).^2,max(rfull(:)).^2),'FontWeight','normal')
       set(findall(gca, 'type', 'text'), 'visible', 'on') 
       fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
       saveCurFigs(get(groot, 'Children'),{'-dsvg'},sprintf('SubspaceNetworks_Indi_%sdim%d',data{1}(1).area_label{i},cur_d),savedir,0); close all
   end  
   
end





end %function 


function [rmat] = SubspaceSimilarityMatrix(data,cur_d)
    
    area_label = data{1}(1).area_label;    
    rmat = NaN(6,14,8,8,8);
    rng('default');
    %get all pairs of subspaces (average across motifs)   
    for cur_rec = 1:6
        fprintf('\nworking on rec %d',cur_rec);
        areas = data{cur_rec}(1).area_label;
        for cur_m = 1:14 
            for f_a = 1:numel(areas)
                for to_a = 1:numel(areas)
                    if to_a>f_a
                        %the following gives you the beta weights for two
                        %different subspaces in all other regions
                        [b,a] = loadBeta(data,cur_d,cur_rec,cur_m,areas(f_a));
                        [bb,aa] = loadBeta(data,cur_d,cur_rec,cur_m,areas(to_a));
                        id = intersect(unique(a),unique(aa));
                        to_idx = find(strcmp(area_label,areas(to_a))==1);
                        f_idx = find(strcmp(area_label,areas(f_a))==1);                        
                        for i = 1:numel(id)
                            temp = LoadActivity(data,cur_rec,cur_m,areas(id(i)));
                            x = temp*b(a==id(i))';
                            y = temp*bb(aa==id(i))';     
                            idx = find(strcmp(area_label,areas(id(i)))==1);
                            rmat(cur_rec,cur_m,f_idx,to_idx,idx) = (corr(x,y,'type','pearson')); %so the connection strength between areas based on how well shared in a third area. 
                        end                        
                    end
                end
            end
        end
    end
    
    


end


function [r,g] = SubspaceSimilarityFull(data,cur_d)
    
    rho = NaN(6,14,28,6);
    grp = cell(6,14,28,6);
    rng('default');
    %get all pairs of subspaces    
    for cur_rec = 1:6
        area_label = data{cur_rec}(1).area_label;
        p = nchoosek(1:numel(area_label),2);
        for cur_m = 1:14 
            for cur_p = 1:size(p,1)
                [b,a] = loadBeta(data,cur_d,cur_rec,cur_m,area_label(p(cur_p,1)));
                [bb,aa] = loadBeta(data,cur_d,cur_rec,cur_m,area_label(p(cur_p,2)));
                id = intersect(unique(a),unique(aa));
                for i = 1:numel(id)
                    temp = LoadActivity(data,cur_rec,cur_m,area_label(id(i)));
                    x = temp*b(a==id(i))';
                    y = temp*bb(aa==id(i))';      
                    rho(cur_rec,cur_m,cur_p,i) = (corr(x,y,'type','pearson'));
                    grp(cur_rec,cur_m,cur_p,i) = area_label(id(i));
                end
            end
        end
    end

    r = rho(:);
    g = grp(:);
    g(isnan(r))=[];
    r(isnan(r))=[];
    r = cellfun(@(x) r(strcmp(g,x)),unique(g),'UniformOutput',0);
    area_label = data{1}(1).area_label;
    %reorder to match area labels (shouldn't actually need this but good to
    %check
    idx = NaN(1,8);
    for i = 1:8
        idx(i) = find(ismember(unique(g),area_label{i}));
    end
    r = r(idx);

end


function [rho,rho_perm,grp,g,n,nn] = SubspaceSimilarity(data,cur_d,targ_area)
rho = NaN(6,14,6);
rho_perm = NaN(6,14,6);
grp = cell(6,14,6);
rng('default');
for cur_rec = 1:6
    for cur_m = 1:14 
        [b,a,n] = loadBeta(data,cur_d,cur_rec,cur_m,targ_area{1});
        [bb,aa,nn] = loadBeta(data,cur_d,cur_rec,cur_m,targ_area{2});
        id = intersect(unique(a),unique(aa));
        [~, area_all] = LoadVariable(data(cur_rec),'rrr_beta','VIS',1); %load betas
        for i = 1:numel(id)
            temp = LoadActivity(data,cur_rec,cur_m,area_all(id(i)));
            x = temp*b(a==id(i))';
            y = temp*bb(aa==id(i))';      
            rho(cur_rec,cur_m,i) = (corr(x,y,'type','pearson'));
            r_perm = NaN(1,1000);
            for j = 1:1000
                r_perm(j) = corr(x(randperm(numel(x),numel(x))),y);
            end
            rho_perm(cur_rec,cur_m,i) = abs(max(r_perm)); %to correct for multiple comparisons
            grp(cur_rec,cur_m,i) = area_all(id(i));
        end
    end
end

%compare all to zero. Also compare all to each other
rho = reshape(rho,size(rho,1)*size(rho,2),size(rho,3));
rho_perm = reshape(rho_perm,size(rho_perm,1)*size(rho_perm,2),size(rho_perm,3));
grp = reshape(grp,size(grp,1)*size(grp,2),size(grp,3));
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load areas
g = NaN(size(grp));
for i = 1:numel(area_all)
   idx = strcmp(grp,area_all{i});
   g(idx)=i;
end
rho = fisherZ(rho);
rho_perm = fisherZ(rho_perm);


end



function ExampleSubspaceGen(data,cur_rec,cur_m,cur_d,source_area,targ_areas)
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load areas
if numel(cur_d)==2
    [b,a,n] = loadBeta(data,cur_d(1),cur_rec,cur_m,targ_areas{1});
    [bb,aa,~] = loadBeta(data,cur_d(2),cur_rec,cur_m,targ_areas{2});
else
    [b,a,n] = loadBeta(data,cur_d,cur_rec,cur_m,targ_areas{1});
    [bb,aa,~] = loadBeta(data,cur_d,cur_rec,cur_m,targ_areas{2});
end
fp = fig_params_cortdynamics;
% temp = LoadActivity(data,cur_rec,cur_m,source_area);
% x = temp*b(a==find(strcmp(area_all,source_area)))';
% y = temp*bb(aa==find(strcmp(area_all,source_area)))';
x = b(a==find(strcmp(area_all,source_area)))';
y = bb(aa==find(strcmp(area_all,source_area)))';


figure; hold on; 
plot(x,y,'marker','.','color',fp.c_area(find(strcmp(area_all,source_area)==1),:),'markersize',fp.markersizebig*1.25,'LineStyle','none');
AddLSline(x,y,x,'k');
rho = corr(x,y);
rng('default')
rho_boot = NaN(1,1000);
for i = 1:1000 
    idx = datasample(1:size(x,1),size(x,1));
    rho_boot(i) = corr(x(idx),y(idx));
end
p = sum(rho_boot<=0)/numel(rho_boot);
ci = [prctile(rho_boot,2.5),prctile(rho_boot,97.5)];
if numel(cur_d)==2
    xlabel(sprintf('%s-%s \n subspace (Dim %d)',source_area,targ_areas{1},cur_d(1)))
    ylabel(sprintf('%s-%s \n subspace (Dim %d)',source_area,targ_areas{2},cur_d(2)))
else
    xlabel(sprintf('%s-%s \n subspace (Dim %d)',source_area,targ_areas{1},cur_d))
    ylabel(sprintf('%s-%s \n subspace (Dim %d)',source_area,targ_areas{2},cur_d))
end
title(sprintf('r %d m %d rho=%0.2f\np=%0.3f ci=%0.3f %0.3f',cur_rec,cur_m,rho,p,ci(1),ci(2)),'FontWeight','normal','fontsize',8)
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 3 4],[2 10 15 15])
end

%plot boxplot
function ContributionPlot(nneu,area_label,data,spreadfactor,xperm)
if nargin <4; spreadfactor=[1,0.75]; end
if nargin <5; xperm=[]; end
fp = fig_params_cortdynamics;
[~, area_name] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
figure; hold on;       
col = fp.c_area(ismember(area_name,area_label),:);
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
%sort
[~,ind] = sort(nanmean(nneu),'descend');
col = col(ind);
CompareViolins(nneu(:,ind)',fp,'plotspread',0,'divfactor',spreadfactor(1),'plotaverage',1,'col',col,'distWidth',spreadfactor(2));
% if ~isempty(xperm)
%    CompareViolins(xperm(:,ind)',fp,'plotspread',0,'divfactor',spreadfactor(1),'plotaverage',0,'col',repmat({[0 0 0]},1,numel(ind)),'distWidth',spreadfactor(2)); 
% end
set(gca,'XTickLabel',area_label(ind),'XTickLabelRotation',45)
fp.FormatAxes(gca);  box on; grid on; 
fp.FigureSizing(gcf,[3 3 6 3],[10 10 14 10])    
ylabel('% of population');

end

%plot boxplot
function MakePlot(b,area,area_label)
fp = fig_params_cortdynamics;
figure; hold on;       
boxplot(b,area,'notch','on')
f = get(get(gca,'children'),'children');
t = get(f,'tag');
idx = find(strcmp(t,'Box')==1);
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', fp.c_area(j,:),'LineWidth',1.5);
end
idx = find(strcmp(t,'Median')==1);    
for j = 1:numel(idx)
    set(f(idx(j)), 'Color', 'k','LineWidth',1.5);
end
set(gca,'XTickLabel',area_label,'XTickLabelRotation',45)
fp.FormatAxes(gca);  box on; grid on; 
fp.FigureSizing(gcf,[3 3 6 3],[10 10 14 10])    
ylabel('|Beta|');
end

function [b,area,area_label] = reorganizeVisBeta(data,cur_d,cur_rec,cur_m,targ_area)
[b, area_label] = LoadVariable(data,'rrr_beta',targ_area,cur_d); %load betas
area_label = area_label(strcmp(area_label,targ_area)==0);
area=LoadVariable(data,'beta_region',targ_area);
%get one rec and one motif
b = abs(squeeze(b(cur_rec,cur_m,:))); %take the absolute values
area = squeeze(area(cur_rec,cur_m,:));
idx = isnan(area);
area(idx)=[];
b(idx)=[];
b = arrayfun(@(n) b(area==n),unique(area(~isnan(area))),'UniformOutput',0);
area = arrayfun(@(n) area(area==n),unique(area(~isnan(area))),'UniformOutput',0);

%reorder by decreasing median beta
[~,reord] = sort(cellfun(@(x) nanmedian(x),b),'descend');
b = b(reord);
area = area(reord);
%renumber area
for i = 1:numel(area)
    area{i} = i*ones(size(area{i},1),1);
end


%condition
b =cat(1,b{:});
area = cat(1,area{:});

end


function [b,area,area_label] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area)
% if cur_d>1 %load the dimension of D1 as a starting point
%     [bdim, ~] = LoadVariable(data,'rrr_beta','VIS',1); %load betas
%     bdim = (squeeze(bdim(cur_rec,cur_m,:)));
% end
[b, area_label] = LoadVariable(data(cur_rec),'rrr_beta',targ_area,cur_d); %load betas
area_label = area_label(strcmp(area_label,targ_area)==0);
area=LoadVariable(data(cur_rec),'beta_region',targ_area);
%get one rec and one motif
b = (squeeze(b(cur_m,:)));
area = squeeze(area(cur_m,:));
idx = isnan(area);
area(idx)=[];
b(idx)=[];
% if cur_d >1
%     bdim(idx)=[];
%     [~,reord] = sort(bdim,'descend');
%     area = area(reord);
%     b = b(reord);
% else
%     %reorder by decreasing median beta
%     [b,reord] = sort(b,'descend');
%     area = area(reord);
% end

%reorder by decreasing median beta
[b,reord] = sort(b,'descend');
area = area(reord);
end


function [b,area,area_label] = loadBeta(data,cur_d,cur_rec,cur_m,targ_area)
[b, area_label] = LoadVariable(data(cur_rec),'rrr_beta',targ_area,cur_d); %load betas
area_label = area_label(strcmp(area_label,targ_area)==0);
area=LoadVariable(data(cur_rec),'beta_region',targ_area);
%get one rec and one motif
b = (squeeze(b(cur_m,:)));
area = squeeze(area(cur_m,:));
idx = isnan(area);
area(idx)=[];
b(idx)=[];
end

function x = LoadActivity(data,cur_rec,cur_m,targ_area)
area_label = data{cur_rec}(cur_m).area_label;
idx = find(strcmp(area_label,targ_area));
x = data{cur_rec}(cur_m).area_val{idx};
%normalize to baseline
x = normalizeToBaseline(x,[1:2],'mean');


%use post stimulus
x = x(:,3:end,:);


%subtract the psth
x = x-nanmean(x,3);


x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
end



% 
% %plot boxplot
% function SimilarityPlot(rho,g,area_label,areas)
% fp = fig_params_cortdynamics;
% figure; hold on;       
% col = fp.c_area(ismember(area_label,areas),:);
% col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
% %sort
% unig = unique(g(~isnan(g)));
% mu = arrayfun(@(n) nanmean(rho(g==n)),unig);
% [~,ind] = sort(mu,'descend');
% col = col(ind);
% unig = unig(ind);
% for i = 1:numel(unig)
%     CompareViolins(rho(g==unig(i))',fp,'plotspread',0,'divfactor',3,...
%         'xpos',i,'plotaverage',1,'col',col(i));
% end
% xlim([0.5 numel(unig)+0.5])
% % figure; hold on;       
% % boxplot((rho),g,'notch','on')
% % f = get(get(gca,'children'),'children');
% % t = get(f,'tag');
% % idx = flipud(find(strcmp(t,'Box')==1));
% % for j = 1:numel(idx)
% %     col = fp.c_area;
% %     col=col(strcmp(area_label,areas{j}),:);
% %     set(f(idx(j)), 'Color', col,'LineWidth',1.5);
% % end
% % idx = find(strcmp(t,'Median')==1);    
% % for j = 1:numel(idx)
% %     set(f(idx(j)), 'Color', 'k','LineWidth',1.5);
% % end
% set(gca,'xtick',1:numel(ind),'XTickLabel',areas(ind),'XTickLabelRotation',45)
% fp.FormatAxes(gca);  box on; grid on; 
% fp.FigureSizing(gcf,[3 3 6 3],[10 10 14 10])    
% ylabel({'Subspace','similarity (rho_z)'});
% 
% end


