function PCANetworks_noProjection(data,ndim)
%Camden - timeless
%load data
if nargin<1; data = LoadSubspaceData('in_grouped'); end
if nargin <2; ndim = 10; end %number of dimensions to use


%full list of areas (recs 1,2,5,6 have all areas; 3 is missing RSP and BFD. 4 if missing RSP)
area_label = data{1}(1).area_label; 

%% Create example figure
cur_rec = 1;
cur_motif = 1; 
% ndim = 10; 
betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,0); 
%row of beta is predictor
area = 1;
fp = fig_params_cortdynamics;
b = cat(2,betas{area,:}); %row 1 is HIPP

%do PCA
[~,~,~,~,pev] = pca(b);

%run a shuffle
rng('default')
pevperm = NaN(numel(pev),1000);
for i = 1:1000
    bperm = b;
    for j = 1:size(b,1)
        temp = b(j,:);        
        bperm(j,:) = temp(randperm(numel(temp),numel(temp)));
    end
    [~,~,~,~,pevperm(:,i)] = pca(bperm);
end

figure; hold on; 
aucperm = NaN(1,1000);
eperm = NaN(1,1000);
for i = 1:1000
   plot(0:numel(pev),cat(1,0,cumsum(pevperm(:,i))),'linewidth',0.1,'color',[0.55 0.55 0.55]); 
   aucperm(i) = trapz([0,cumsum(pevperm(:,i))']/numel([0,pevperm(:,i)']));
   eperm(i) = ED(pevperm(:,i));
end
plot(0:numel(pev), cat(1,0,cumsum(pev)),'color','k','linewidth',1.5);
ylim([0 100]);
ylabel('% Explained Variance');
xlabel({'Dimension','(7 subspace x 10 dimensions)'});
set(gca,'xtick',[0,35,70],'ytick',[0:25:100])
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 2.5 3.25],[10 10 12 10])

% get global statistics
e = ED(pev);
auc = trapz([0,cumsum(pev)']/numel([0,pev']));

figure; hold on; 
histogram(eperm,'BinWidth',1,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.8)
yvals = get(gca,'ylim');
plot([e,e],yvals,'color','k','linewidth',2,'linestyle',':');
fp.FormatAxes(gca);  box on; grid on
xlim([5 20])
ylim(yvals);
fp.FigureSizing(gcf,[3 3 1.75 1.75],[10 10 12 10])
p = sum([eperm,e]<=e)/numel([eperm,e]);

title(sprintf('auc %0.3f%% ed=%0.3f \n auc %0.3f%% ed=%0.3f \n p=%0.3f',auc,e,nanmean(aucperm),nanmean(eperm),p),'FontWeight','normal')


%% Get statistic across all fits
if exist('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\PCANetworks\PCA_data_10dim_EffectiveDim.mat','file')
    load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\PCANetworks\PCA_data_10dim_EffectiveDim.mat');
else %recompile
    %%
    nrec = numel(data); %num  recordings
    nmotif = size(data{1},2); %num motifs
    stats = cell(nrec,nmotif); %preallocate
    for cur_rec = 1:nrec %loop through recordings
        fprintf('\nworking on rec %d',cur_rec);
        for cur_motif = 1:nmotif %loop through motifs    
            %load betas            
             betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,0);             
             stats{cur_rec,cur_motif} = getSummaryStatsPerArea(betas,1);
        end %motif loop    
    end %rec loop  

    stats = cat(1,stats{:});
    stats = cat(1,stats(:));
    idx = arrayfun(@(n) stats(n).area,1:size(stats,1),'UniformOutput',0);
    idx = cellfun(@ isempty, idx);
    stats(idx)=[];
    
    save('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\PCANetworks\PCA_data_10dim_EffectiveDim.mat');
end
%remove empties

%% Make summative plots
fp = fig_params_cortdynamics;

%plot the delta AUC/Raw across areas
d = NaN(1000,8);
for i = 1:8
   idx = find(cat(1,stats(:).area)==i);
   temp = cat(1,stats(idx).e);
   tempperm = nanmean(cat(1,stats(idx).eperm),2);
   d(:,i) = pairedBootstrap(temp./tempperm,@nanmean);   
end
[~,idx] = sort(nanmean(d),'descend');
d = d(:,idx)*100;
area_name = area_label(idx);

%reorder
col = fp.c_area; 
col = col(idx,:); 
col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);

figure; hold on; 
CompareViolins(d',fp,'plotspread',0,'divfactor',0.5,'plotaverage',1,'col',col,'distWidth',0.75);
set(gca,'XTickLabel',area_name,'XTickLabelRotation',45)
% yval = get(gca,'ylim');
set(gca,'ylim',[35 60]);
% set(gca,'ylim',[floor(yval(1)),ceil(yval(2))]);
ylabel({'Sharing D_{true}/D_{shuf}'});
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 7 4],[10 10 12 10])

%%get global statistics
temp = cat(1,stats(:).e);
tempperm = nanmean(cat(1,stats(:).eperm),2); 
[~,bootstat] = pairedBootstrap((temp./tempperm)*100,@nanmean);
title(sprintf('mu %0.3f%% CI: %0.3f - \n%0.3f p=%0.3f',bootstat.mean,bootstat.ci(1),bootstat.ci(2),bootstat.p),'FontWeight','normal')

%do significance testing: 
pval = NaN(8,8);        
figure; hold on;        
for i = 1:size(pval,1)
    for j= 1:size(pval,2)
       temp = d(:,i)-d(:,j);
       pval(i,j) = sum(temp>0)/numel(temp(~isnan(temp)));
       plot(i,j,'marker','.');
       text(i,j,sprintf('%0.3f',pval(i,j)),'FontWeight','normal','fontsize',10);
    end
end
set(gca,'xtick',[1:8],'XTickLabel',area_name,'ytick',[1:8],'yticklabel',area_name);

%get the mean and CI for each
mean(d,1)
prctile(d,2.5,1)
prctile(d,97.5,1)


%% Get all correlations across all beta values
p = cell(6,14); %preallocate
for cur_rec = 1:6 %loop through recordings
    fprintf('\nworking on rec %d',cur_rec);
    for cur_motif = 1:14 %loop through motifs    
        %load betas            
        betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,0);             
        p{cur_rec,cur_motif} = SigRelationship(betas);
    end %motif loop    
end %rec loop




end %function end


function p = SigRelationship(betas)
 p = NaN(8,70,70);
 for cur_a = 1:8
    b = cat(2,betas{cur_a,:}); 
    for i = 1:size(b,2)
        for j = 1:size(b,2)
            if j>i
                rho = corr(b(:,i),b(:,j));
                if rho>0 %since direction is arbitrary
                    [~,p(cur_a,i,j)] = corr(b(:,i),b(:,j),'tail','right');
                else
                    [~,p(cur_a,i,j)] = corr(b(:,i),b(:,j),'tail','left');
                end
%                 rboot = NaN(1,1000);
%                 rng('default');
%                 for q = 1:1000 
%                     idx = datasample(1:size(b(:,i),1),size(b(:,i),1));
%                     rboot(q) = corr(b(idx,i),b(idx,j));                            
%                 end
%                 if rho>0 %since direction is arbitrary
%                     p(cur_a,i,j) = sum(rboot<0)/numel(rboot);
%                 else
%                     p(cur_a,i,j) = sum(rboot>0)/numel(rboot);
%                 end
            end
        end
    end
 end
 p = p(~isnan(p));
end
%%
function stats = getSummaryStatsPerArea(betas,doperm)
    stats = [];
    for i = 1:size(betas,1)
        %loop through each row (predictor area)
        b = cat(2,betas{i,:});
        if ~isempty(b)
            stats(i).area = i;
            [~,~,~,~,pev] = pca(b);
            stats(i).e = ED(pev);
            stats(i).auc = trapz([0,cumsum(pev)']/numel([0,pev']));
            
            if doperm                
                pevperm = NaN(numel(pev),1000);
                for z = 1:1000
                    bperm = b;
                    for j = 1:size(b,1)
                        temp = b(j,:);        
                        bperm(j,:) = temp(randperm(numel(temp),numel(temp)));
                    end
                    [~,~,~,~,pevperm(:,z)] = pca(bperm);
                end                
                aucperm = NaN(1,1000);
                eperm = NaN(1,1000);
                for z = 1:1000
                   aucperm(z) = trapz([0,cumsum(pevperm(:,z))']/numel([0,pevperm(:,i)']));
                   eperm(z) = ED(pevperm(:,z));
                end     
                stats(i).aucperm = aucperm;
                stats(i).eperm = eperm;
                stats(i).pval = sum([eperm,stats(i).e]<=stats(i).e)/numel([eperm,stats(i).e]);
            end
        end        
    end

end


function [betas,targidx] = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj)
    %retruns betas which is an 8x8 cell array. Each cell cointains matrix
    %of betas weights of a region (nuerons x dimensions)
    % useProj, then returns (timepoints x trials) x dimensions matrix of
    % activtiy in a region projected along those beta weights
    %rows of beta are predictor regions. columns are targets
    %i.e., row 1 of beta column 2 is (neurons x dimensions) betas for Hipp predicting area MOs
    %Gutchecks: 
    %betas should have nothing on the diagnol
    %betas for rec three rows and columns 4/6 should be empty and rec 4 row/column 3

    b = data{cur_rec}(cur_motif).rrr_B; %each cell is beta weights of n-1 regions predicting each target region
    cur_areas = data{cur_rec}(cur_motif).area_label; %the list of target regions
    
    %just keep desired dimensions
    b = cellfun(@(x) x(:,1:ndim),b,'UniformOutput',0);     
    
    %reorganize betas so each cell is the betas for a given area when predicting all other areas
    g = data{cur_rec}(cur_motif).grouping; %each cell is region of beta weights of n-1 regions predicting each target region             
    betas = cell(numel(area_label),numel(area_label));
    targidx = ismember(area_label,cur_areas); %get true index for targets (out of 8 areas, to match across recs)
    for cur_a = 1:numel(cur_areas) 
        idx = strcmp(area_label,cur_areas(cur_a)); %get true index (out of 8 areas, to match across recs)        
        if ~isempty(idx)
            betas(idx,targidx) = arrayfun(@(n) b{n}(g{n}==cur_a,:),1:numel(g),'UniformOutput',0);  
        end
    end  
    
    %if projecting 
    if useProj == 1
       y = data{cur_rec}(cur_motif).area_val;
       y = cellfun(@(x) normalizeToBaseline(x,[1:2],'mean'), y,'UniformOutput',0);%normalize to baseline
       y = cellfun(@(x) x(:,3:end,:), y,'UniformOutput',0);%use non-baseline
       y = cellfun(@(x) x-nanmean(x,3), y,'UniformOutput',0);%subtract the psth
       y = cellfun(@(x) reshape(x,[size(x,1),size(x,2)*size(x,3)])', y,'UniformOutput',0);%reshape into trial/timepoint x neuron   
       
       %reshape to match predictors in beta
       yy = cell(8,1);
       yy(targidx)=y;
       yy = repmat(yy,1,8);
       %temporily replace emptys with nan for multiplication
       badidx = cellfun(@isempty,betas);
       yy(badidx) = {NaN};
       betas(badidx) = {NaN};
       betas = cellfun(@(x,y) x*y, yy,betas,'UniformOutput',0);
       betas(badidx) = {[]};
       
    else
        badidx = cellfun(@isempty,betas);        
        betas(badidx) = {[]};
    end
end


function e = ED(pev)
if sum(pev)>1
    pev=pev/100;
end
    
e = 1/sum(pev.^2);

end



% 
% function [w,b] = withinbetweenEnt(coef,ndim)
%     get_ent = @(p) -sum(p.*log2(p));
%     
%     w = NaN(1,7*ndim);
%     b = NaN(1,7*ndim);
%     for j = 1:size(coef,2)
%         p = abs(coef(:,j));
%         %get the shape of the distribution across PCA dimensions
%         byarea = reshape(p,size(p,1)/ndim,ndim);
%         narea = size(byarea,1);
%         if narea>ndim
%             %entropy within target areas
%             within = nanmean(arrayfun(@(n) get_ent(byarea(n,:)),1:size(byarea,1)));
%             %subsample across areas
%             between = nanmean(arrayfun(@(n) get_ent(byarea(randperm(narea,ndim),n)),1:size(byarea,2)));
%         elseif narea<ndim
%             %entropy within target areas
%             within = nanmean(arrayfun(@(n) get_ent(byarea(n,randperm(ndim,narea))),1:size(byarea,1)));
%             %subsample across areas
%             between = nanmean(arrayfun(@(n) get_ent(byarea(:,n)),1:size(byarea,2)));
%         else %
%             %entropy within target areas
%             within = nanmean(arrayfun(@(n) get_ent(byarea(n,:)),1:size(byarea,1)));
%             %subsample across areas
%             between = nanmean(arrayfun(@(n) get_ent(byarea(:,n)),1:size(byarea,2)));
%         end
%         w(j) = within;
%         b(j) = between;
%     end
%     
% 
% end
% 
% function [dauc] = nullAUC(b)
%     %null 
%     dauc = NaN(1,1000);
%     for i = 1:1000
%         bnull = b;
%         for j = 1:size(b,2)
%             idx = reshape(1:size(b,1),12,size(b,1)/12);
%             idx = idx(:,randperm(size(idx,2),size(idx,2)));
%             bnull(:,j) = bnull(idx(:),j);
% %             bnull(:,j) = bnull(randperm(size(b,1),size(b,1)),j); %full random            
%         end
%         [dauc(i),~,~] = getAUC(bnull,0);
%     end
% end
% 
% 
% function [dauc,auc,auc_raw,coef] = getAUC(b,oldflag)  
%     if nargin <2; oldflag=1; end
%     x = 100*sort(nanvar(b)/sum(nanvar(b)),'descend'); %raw variance
%     [coef,~,~,~,y] = pca(b); %pca variance    
%     if oldflag   
%         % do it the old way (subtraction)
%         %compare AUC
%         x = [0,cumsum(x)];
%         y = [0,cumsum(y')];                
%         auc_raw = trapz(x/numel(x));
%         auc = trapz(y/numel(y));
%         dauc = auc-auc_raw;  
%     else
% 
%         %do the new way (compare)
%         a=trapz([0,cumsum(x)], [0,cumsum(y)'])/100;
% 
%         auc_raw =[];
%         auc=[];
%         dauc=a;
%     end
% end










