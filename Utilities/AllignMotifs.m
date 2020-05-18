function X_alligned = AllignMotifs(X,core_comm_idx,lags,clust_id,lag_mat)

[P, N, Z] = size(X);

%pad X
X_pad = NaN(P,N,Z*3);
for i =1:N
    X_pad(:,i,:) = cat(2, zeros(P,Z), squeeze(X(:,i,:)), zeros(P,Z));     %Zero pad
end
X_alligned = X_pad;
for i = 1:numel(unique(clust_id))
    template_idx = core_comm_idx{i}(randperm(numel(core_comm_idx{i}),1)); %get a template pattern from that cluster; 
    cluster_idx = find(clust_id==i);    
    
    %allign all to the template use the lags
    for j = 1:numel(cluster_idx)
        shift_val=lag_mat(template_idx,cluster_idx(j));        
        X_alligned(:,cluster_idx(j),:) = circshift(squeeze(X_pad(:,cluster_idx(j),:)),lags(shift_val),2);
    end
end
    
end