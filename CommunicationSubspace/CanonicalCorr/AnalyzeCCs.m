function AnalyzeCCs(folder,rec_name)
%Camden - timeless
%@input folder is the folder containing the analyzed CCA

%grab data
rec_name = '';
[fn,~] = GrabFiles([rec_name '\w*RRR\w*.mat'],0,{folder}); 
data = cellfun(@(x) load(x,'motif','a','b','U','V','r','paired_areas','area_label','best_idx','t','aTheta_xv','bTheta_xv','pval'),fn);

%% main






% %get the strength
% [cv_map,x] = cvStrengthMap(data,'r_first');
% for i = 
% 
% 


%
%for each subspace, grab the strongest motif and compare to the 

%Area list
area_label = data(1).area_label;
paired_areas = data(1).paired_areas;
[cort_ss,cort_idx] = GetCorticalSubspaces(paired_areas, area_label);

a = [data.a]; 
idx = cellfun(@(x) isempty(x),a);
a(idx) = {NaN};
a = cellfun(@(x) x(:,1),a,'UniformOutput',0);
%loop through each subspace and compare the angle between coeficients
for i = 1:size(cort_idx,1)   
   temp = a(cort_idx(i),:);   
   %remove empty ones
   temp(cellfun(@(x) sum(isnan(x))==1,temp,'UniformOutput',1))=[];
   temp = [temp{:}];
   cstat = SubspaceConsistency(temp,temp,'angle',1);
   rng('default')
   cstat_perm = NaN(1,1000);
   for j = 1:1000
      perm_idx = arrayfun(@(n) randperm(size(temp,1),size(temp,1)),1:size(temp,2),'UniformOutput',0);
      perm_idx = cat(1,perm_idx{:})';
      cstat_perm(j) = nanmean(SubspaceConsistency(temp(perm_idx),temp,'angle',1));
   end
end

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


%note camden, this is doing the cca per timepoint
%get the similarity in angle between all significant timepoints
%(variability suggests not a valid cc)
%show a heatmap of when the significant cvs occur
idx = data(1).t{1};
for i = 1:numel(data)  
    %get strength
    r = data(i).r;
    temp_idx = cellfun(@(x) isempty(x),r);
    r(temp_idx) = [];  
    r = cellfun(@(x) x(1),r,'UniformOutput',1);  
    r = repmat(r,1,2)';
    r = r(:);

    %get number per time
    temp = data(i).best_idx;
    temp_idx = cellfun(@(x) isempty(x),temp);
    temp(temp_idx) = [];  
    temp = cellfun(@(x) x(1),temp,'UniformOutput',1);
    
    temp = arrayfun(@(x) idx(:,x),temp,'UniformOutput',0);
    
    temp = cat(1,temp{:});
    
    figure; 
    plot(temp,r,'linestyle','none','marker','x')
    title(sprintf('%d',i));
end

%compare the strength of the best CC for each motif 
[cv_map,x] = cvStrengthMap(data,'r_first');
fh = visualizeStrengthMap(cv_map,data);

cv_map = cvStrengthMap(data,'r_all_sum');
fh = visualizeStrengthMap(cv_map,data,[]);

%similarity between upper and lower triangles


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