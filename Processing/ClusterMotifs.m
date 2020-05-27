function [basis_motifs, kval, ovr_q, idx_louvain, idx_knn, tcorr_mat, handles] = ClusterMotifs(file_list,gp)
%camden macdowell - timeless
%file_list is the full path of each file with motifs to cluster. motifs are
%contained in variable 'w' as a pxl x motif x time tensor

if nargin <2; gp = general_params; end


%load motifs
fprintf('\n\t Loading Motifs from Individual Chunks...')
motifs = cell(1,numel(file_list));
nanpxs = cell(1,numel(file_list));
for cur_file = 1:numel(file_list)   
    if mod(cur_file,floor(numel(file_list)*.25))==0
        fprintf('\n%d %% Complete',round(cur_file/numel(file_list)*100));
    end
    temp = load(file_list{cur_file},'w','nanpxs');
    nanpxs{cur_file} = temp.nanpxs;
    motifs{cur_file} = temp.w(:,nanvar(squeeze(sum(temp.w,1)),[],2)>0,:); %only use filled motifs
end

%Use only pixels that are not masked for all motifs (e.g if clustering across different animals)
if numel(unique(cellfun(@numel, nanpxs,'UniformOutput',1)))>1 %if different nan patterns
    mask_stack = ones(size(motifs{1},1)+numel(nanpxs{1}),numel(nanpxs));
    for cur_motif = 1:numel(nanpxs) %reconstruct full mask for each    
       temp = conditionDffMat(ones(size(motifs{1},1),1)',nanpxs{cur_motif},[],[gp.pixel_dim,1]);
       mask_stack(:,cur_motif) = temp(:);       
    end
    %get the global nanmask
    nanpxs = find(any(isnan(mask_stack),2)==1)';    
else
    nanpxs = nanpxs{1};
end
  
%concatenate
motifs = cat(2,motifs{:});
[P, N, Z] = size(motifs);

%optional 3D gaussian smooth
if ~isempty(gp.m_smooth_kernel)
    fprintf('\n\t Smoothing Motifs....');
    for cur_motif = 1:N
        temp = conditionDffMat(squeeze(motifs(:,cur_motif,:))',nanpxs,[],[gp.pixel_dim,size(motifs,3)]);
        bad_pxls = isnan(temp); %get a mask of the nan pixels so that image size doesn't change with smoothing
        temp(bad_pxls)=0;
        temp = imgaussfilt3(temp,gp.m_smooth_kernel);
        temp(bad_pxls)=NaN;
        temp = conditionDffMat(temp)';
        motifs(:,cur_motif,:) = reshape(temp,P,1,Z);
    end 
end

%could include an optional renormalization

%compute maximum temporal xcorrelation between motifs
[tcorr_mat, lag_mat, lags] = TemporalXCorrTensor(motifs,gp.m_maxshift);

%Fit Phenograph Number of Neighbors
[kval, ~] = FitPhenoK(tcorr_mat,'k_range',gp.pheno_k_range,'louvain_restarts',gp.pheno_louvain_restarts,'genfigs',1,'number_resamples',gp.pheno_num_resamples);

%Cluster
[idx_louvain, idx_knn, ovr_q] = PhenoCluster(tcorr_mat,'k',kval,'louvain_restarts',gp.pheno_louvain_restarts,'Verbose',0);   
if size(idx_louvain,2) > 1; idx_louvain = idx_louvain(:,end); end %just take the initital clustering. Phenograph tends to really overcluster at the higher levels. 

%visualize cross correlation matrix
Plot_OrderedSimilarityMatrix(tcorr_mat,idx_louvain);

%Get core community to average for basis motifs
[core_comm_idx, ~] = CoreCommunity(idx_louvain,idx_knn,gp.m_community_fraction); 

%Allign motifs in each cluster to one of the core community members 
motifs_alligned = AllignMotifs(motifs,core_comm_idx,lags,idx_louvain,lag_mat);

%center shift
motifs_alligned = helper.shiftW(motifs_alligned);

%compute basis motifs
basis_motifs = NaN(P,numel(core_comm_idx),size(motifs_alligned,3));
for i = 1:numel(core_comm_idx)   
    basis_motifs(:,i,:) = nanmean(motifs_alligned(:,core_comm_idx{i},:),2);
end

%optional removal of the padded pixels
if gp.m_removepad
   basis_motifs = basis_motifs(:,:,Z+1:Z*2);
else %just remove totally empty regions
   basis_motifs = basis_motifs(:,:,nanvar(squeeze(sum(basis_motifs,1)),[],1)>eps);
end

handles = get(groot, 'Children');

end


















        