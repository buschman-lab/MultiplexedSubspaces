function AnalyzeCCs(folder,rec_name)
%Camden - timeless
%@input folder is the folder containing the analyzed CCA

%grab data
[fn,~] = GrabFiles([rec_name '\w*fullrec\w*.mat'],0,{folder}); 
data = cellfun(@(x) load(x,'motif','a','b','U','V','r','r_norm','paired_areas','area_label'),fn(1));
motif_cvs = cellfun(@(x) load(x,'motif_cvs'),fn);

%% main





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