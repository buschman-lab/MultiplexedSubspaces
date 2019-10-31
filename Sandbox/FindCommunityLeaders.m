%find the motifs with the most nearest neighbors in the same graph

rel = Idx_knn;
for i = 1:size(Idx_knn,1)
    for j = 1:size(Idx_knn,2)
        rel(i,j) = idx_louvain(Idx_knn(i,j),2);
    end
end

%loop through each row and find those that best correspond to their
%community

numcom = NaN(size(rel,1),1);
for i = 1:size(rel,1)
    numcom(i) = sum(rel(i,:)==rel(i,1));
end

%Get only the ones that are most within community
val = unique(idx_louvain(:,2));
best_idx = [];
for i = 1:numel(val)
    cur_com_id = find(idx_louvain(:,2)==val(i));
    [cur_com, loc] = sort(numcom(cur_com_id),'ascend');   
    %get the top 10% in each community
    x = floor(0.10*length(cur_com));
    loc = loc(end-x:end);
    best_idx{i} = cur_com_id(loc);
    com_leadership{i} = cur_com(end-x:end);
end

%Get the alligned version and get the average for each community; 
basis = zeros(4624,numel(val),29);
for i = 1:numel(val)
    temp = W_smooth_alligned(:,best_idx{i},:);
    basis(:,i,:) = nanmean(temp,2);
    temp = squeeze(nanmean(temp,2));
    temp = reshape(temp,[68 68 29]);
    figure; 
    for j = 1:29
        title(sprintf('motif %d',i));
        cla; imagesc(temp(:,:,j),[0 0.5]); pause(); 
    end
    close; 
end
  
 
    
    
    
    
    
    
    