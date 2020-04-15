function VarianceOfMotifsCapturedByBasis()
%Camden MacDowell - timeless

%load data
idx = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\ClusteredDPs_Paramset.mat','idx_louvain','Idx_knn');
idx_louvain = idx.idx_louvain;
Idx_knn= idx.Idx_knn;
W_smooth_alligned = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\AverageDPs_suppliment_1','W_smooth_alligned');
W_smooth_alligned = motifs.W_smooth_alligned;
%Unsmoothed



%First do smoothed
[core_community_indices,core_community_values] = CoreCommunity(idx_louvain(:,end),Idx_knn,0.10);
%Get the alligned version and get the average for each community; 
W_clust_smooth = zeros(size(W_smooth_alligned,1),numel(core_community_indices),size(W_smooth_alligned,3));
for i = 1:numel(core_community_indices)   
    W_clust_smooth(:,i,:) = nanmean(W_smooth_alligned(:,core_community_indices{i},:),2);
end

%compute the explained variance for each group by it's basis
val = unique(idx_louvain);
expvar_all = cell(1,numel(val));
for i = 1:numel(val)
    cur_com_id = find(idx_louvain(:,end)==val(i));
    basis = squeeze(W_clust_smooth(:,i,:));    
    temp_expvar = NaN(1,numel(cur_com_id));
    for cur_motif = 1:numel(cur_com_id)
        motif = squeeze(W_smooth_alligned(:,cur_com_id(cur_motif),:));
        resid = motif-basis;
        temp_expvar(cur_motif) = CalculateExplainedVariance(motif,resid);
    end
    expvar_all{i} = temp_expvar;    
end

expvar_median = cellfun(@nanmedian,expvar_all);
expvar_ci_median = cellfun(@(x) bootci(1000,@nanmedian,x),expvar_all,'UniformOutput',0);
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ManuscriptRevisionFigures_currentbio';
save([savedir filesep 'motiffitwithbasis.mat'],'expvar_median','expvar_all','expvar_ci_median');

%save off

%
