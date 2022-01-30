function RemoveNoiseClusters()
%Basis motif figure/within animal clustering
%Camden MacDowell - timeless

cd('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscovery_temp');
basis = 'Mouse_basis_motifs.mat';
mouse_fn = {'Mousepermouse_basis_motifs331.mat','Mousepermouse_basis_motifs332.mat','Mousepermouse_basis_motifs334.mat'};

load(basis); 
unique_cluster_idx = unique(cluster_idx);

%remove motifs that are not shared across all animals
mouse_clusters = load(mouse_fn{1},'cluster_idx_all');
mouse_clusters = mouse_clusters.cluster_idx_all; 
not_shared = cellfun(@(x) find(ismember(1:numel(unique_cluster_idx),unique(x))==0),mouse_clusters,'UniformOutput',0);
not_shared = unique([not_shared{:}]);

%get balance in each cluster
clust_contrib = NaN(numel(mouse_fn),numel(unique_cluster_idx));
for i = 1:numel(unique_cluster_idx)
    clust_contrib(:,i) = cellfun(@(x) sum(x==unique_cluster_idx(i)),mouse_clusters,'UniformOutput',1);
end

shared = find(~ismember(1:numel(unique_cluster_idx),not_shared));

%contributed to less than x number of epochs 
bad_motifs = ismember(1:numel(unique_cluster_idx),not_shared);
bad_motifs = bad_motifs + sum(clust_contrib)<25;

% max(clust_contrib)./sum(clust_contrib)
arrayfun(@(n) fprintf('\nCluster %d | # motifs = %d',n,sum(clust_contrib(:,n))),shared,'UniformOutput',1)

fprintf('\n\tNumber of nonshared motifs = %d from %d clusters\n',sum(clust_contrib(:,not_shared),'all'),numel(not_shared))

%manual removal 
w_good = W_basis; 
w_good(:,bad_motifs,:)=[];

close all
for i = 1:size(w_good,2)
    temp = reshape(w_good(:,i,:),[68 68 size(w_good,3)]);
    for j = 1:size(temp,3)
        imagesc(temp(:,:,j),[0 prctile(temp(:),99.99)]);
        title(sprintf('%d %d',i,j));
        pause();
    end
end

temp = input('list of obvious noise clusters:    ');

bad_motifs(temp)=true;

fprintf('\n\tNumber of nonshared motifs = %d from %d clusters\n',sum(clust_contrib(:,bad_motifs),'all'),sum(bad_motifs))

noise_clusters = find(bad_motifs==1);
save(basis,'noise_clusters','-append');
for i = 1:numel(mouse_fn)
   save(mouse_fn{i},'noise_clusters','-append');
end

end %function 


