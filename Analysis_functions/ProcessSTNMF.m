function stnmf_out = ProcessSTNMF(stnmf,original_data)
%Camden MacDowell - timeless
%this takes data structure stnmf and computes expvar of the fit as a
%function of ordered dimensions in A (ordered by decreasing weighting). 
    
max_dim = 50;
pev = NaN(numel(stnmf),max_dim);
for cur_rec = 1:numel(stnmf)    
    fprintf('\n\t Working on rec %d of %d',cur_rec,numel(stnmf));
    x = stnmf(cur_rec).fit_full;  %get full fit     
  
    [coef, score, ~, ~, ~, mu] = pca(x.A'); %pca is columnwise so this is spatial modes of A.     
    %Reconstruct the sparse A
    for cur_dim = 1:max_dim        
        sparse_a = repmat(mu,size(score,1),1) + score(:,1:cur_dim)*coef(:,1:cur_dim)';   
        sparse_a = sparse_a';
        X_recon = stNMF_Reconstruct(sparse_a,x.Wt,x.Ws);
        X_recon = X_recon';
%         pev(cur_rec,cur_dim) = CalculateExplainedVariance(original_data{cur_rec},...
%             original_data{cur_rec}-X_recon)*100;
        pev(cur_rec,cur_dim) = CalculateExplainedVariance(fullx,...
            fullx-X_recon)*100;
    end
    
end
    
    
%     sparse_a = zeros(size(x.A));
%     temp_a = x.A; %used to iteratively zero out maxes
%     
%     for cur_dim = 1:max_dim        
%         [coord(cur_dim,3), idx] = max(temp_a(:)); %get index
%         [coord(cur_dim,1),coord(cur_dim,2)] = ind2sub(size(temp_a),idx); %convert to coordinates
%         temp_a(coord(cur_dim,1),coord(cur_dim,2))=0; %zero out for next iteration                
%         sparse_a(coord(cur_dim,1),coord(cur_dim,2)) = coord(cur_dim,3); %populate with top relationship
%         %compute the explained variance
%         X_recon = stNMF_Reconstruct(sparse_a,x.Wt,x.Ws);
%         X_recon = X_recon';
%         pev(cur_rec,cur_dim) = CalculateExplainedVariance(X,...
%             X-X_recon)*100;    
%     end %dimension loop 