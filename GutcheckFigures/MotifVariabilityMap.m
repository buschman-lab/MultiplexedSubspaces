function MotifVariabilityMap(file_path)
% file_path = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery';
% basis_name = 'Mouse_basis_motifs.mat'

file_list = GrabFiles(['Mouse', '\w*chunk\w*.mat'],0,{file_path});
%load all the basis motifs
temp = cellfun(@(x) load(x,'w','nanpxs'),file_list,'UniformOutput',0);
W_orig = cellfun(@(x) x.w, temp,'UniformOutput',0);
nanpxs_all = cellfun(@(x) x.nanpxs,temp,'UniformOutput',0);

%load the clustering information
file_list = GrabFiles(basis_name,0,{file_path});
load(file_list{1},'noise_idx','lag_mat','cluster_idx','core_comm_idx','lags','cluster_idx_all','idx_knn','shift')

%create basis motifs per mouse
mouseid = MouseNumFromPath(file_list,'Mouse_'); unique_mice = unique(mouseid);
parameter_class = 'general_params_corticaldynamics';
gp = loadobj(feval(parameter_class));
gp.clust_community_fraction = 1; gp.clust_removepad=1;
opts = gp;
unique_mice = unique(mouseid);

%match the cluster_idx
mouseid_all = arrayfun(@(n) ones(size(W_orig{n},2),1)*mouseid(n),1:numel(mouseid),'UniformOutput',0);
mouseid_all = cat(1,mouseid_all{:});


nanpxs = nanpxs_all;
%loop through mice and combine the motifs for that animal
W_mouse = cell(1,numel(unique_mice));
core_comm_size = cell(1,numel(unique_mice));
cluster_idx_all = cell(1,numel(unique_mice));
for cur_mouse = 1:numel(unique_mice)
    fprintf('\n\t Working on mouse %d of %d',cur_mouse,numel(unique_mice));
    %confirm that all nanpxs masks are the same for this animal (required)
    try temp_nanpxs = cat(1,nanpxs{mouseid==unique_mice(cur_mouse)}); catch; error('nanpxs masks for the current animal are not consistent across recordings'); end
    temp_nanpxs = temp_nanpxs(1,:);
    
    if iscell(W_orig) %combine mutliple CNMF Fits stored in a cell array
        W = cat(2,W_orig{mouseid==unique_mice(cur_mouse)});
    end

    %Remove empty motifs (legacy compatibility)
    [W,~,indempty] = RemoveEmptyMotifs(W);    
    if ~isempty(indempty)
        idx = find(mouseid_all==unique_mice(cur_mouse)); %indices of cur mouse
        mouseid_all(idx(indempty))=[];
    end
    
    [~, ~, L] = size(W);
    
    %get the clustering info for this animal 
    idx = find(mouseid_all==unique_mice(cur_mouse)); %indices of cur mouse
    temp_cluster_idx = cluster_idx(idx);
    temp_idx_knn = idx_knn(idx,:);
    temp_idx_knn(~ismember(temp_idx_knn,idx))=NaN; %only keeps neighbors from same animal 
    temp_idx_knn = temp_idx_knn-min(temp_idx_knn)+1; %adjust for the shift of a middle animal
    temp_lag_mat = lag_mat(idx,idx);   
    temp_shift = shift(idx);
    
    %optional 2D or 3D gaussian smooth. Reccomended for noisy data. 
    if ~isempty(opts.clust_smooth_kernel)
        W_smooth = GaussianSmoothTensor(W,opts.clust_smooth_kernel,opts.originaldimensions,temp_nanpxs,opts.clust_nobleed);
    else
        W_smooth = W; 
    end
        
    %Reconstruct full W. This is used for averaging to make the basis motifs. 
    W = MaskTensor(W,temp_nanpxs,[opts.originaldimensions(1)*opts.originaldimensions(2),size(W,2),size(W,3)]); 

%     fprintf('\n\tGenerating Basis Motifs for mouse %d',unique_mice(cur_mouse))
%     %Get core community to average for basis motifs
%     if numel(opts.clust_community_fraction)>1 %find the best core_community_fraction
%         fprintf('\n\t Autofitting Community Fractions');
%         [core_comm_idx, core_comm_size{cur_mouse}] = AutoFitCommunityFraction(temp_cluster_idx,temp_idx_knn,opts,W_smooth,temp_lag_mat,lags);               
%     else %just take the set value
%         fprintf('\n\t Using Set Community Fraction');
%         [core_comm_idx, core_comm_size{cur_mouse}] = CoreCommunity(temp_cluster_idx,temp_idx_knn,opts.clust_community_fraction); 
%     end

    %Allign all motifs to the index of the original recording
    fprintf('\n\tAlligning W')
%     W_alligned = AllignW(W,core_comm_idx_all,lags,cluster_idx,lag_mat);
%     %alligned within mice
    W_alligned = AllignW_withinMice(W,lags,temp_shift);

    %compute the variability per pixel in each basis motif
    for i = 1:numel(unique(temp_cluster_idx))
        var_map = [];
        for j = 1:size(W_alligned,3)
           temp = W_alligned(:,temp_cluster_idx==i,j);
           resid = temp-nanmean(temp,2);
           %how large resid relative to mean
           var_map(:,j) = nanmedian(abs(resid)./nanmean(temp,2),2);
       end
    end
       
       
       
       nanmean(W_alligned(:,core_comm_idx{i},:),2); 
    end
    
    %compute basis motifs using just the alligned motifs for this animal
    fprintf('\n\tComputing Basis Motifs')
    W_basis = NaN(size(W,1),numel(core_comm_idx),size(W_alligned,3));
    for i = 1:numel(core_comm_idx)   
        W_basis(:,i,:) = nanmean(W_alligned(:,core_comm_idx{i},:),2);
    end

    %optional removal of the padded pixels
    if opts.clust_removepad
       W_basis = ShiftW(W_basis); %center shift before optional 
       W_basis = W_basis(:,:,L+1:L*2);
    else %just remove totally empty regions
       W_basis = W_basis(:,:,nanvar(squeeze(sum(W_basis,1)),[],1)>eps);
    end   
    W_mouse{cur_mouse} = W_basis;
    
    %save off the temp clustering index for future saving
    cluster_idx_all{cur_mouse} = temp_cluster_idx;
end


