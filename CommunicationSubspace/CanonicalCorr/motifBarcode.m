function motifBarcode(data)
%Camden - timeless
%takes structure data, gets the subspace barcode for each motif (can be
%computed different ways) and generates motif x motif similarity matrix


%just get cortical areas
area_label = data(1).area_label;
paired_areas = data(1).paired_areas;
[~,cort_idx] = GetCorticalSubspaces(paired_areas, area_label);

[~,x] = cvStrengthMap(data,'r_first');
x(:,2)=[];
x = x(cort_idx,:);
x(isnan(x))=0;
figure; hold on; 
plot(x)
ylabel('subspace stength')
xlabel('subspace'); 
title('subspace strength colored motif (no normalization')

x = x-nanmean(x,2);
figure; hold on; 
plot(x)
ylabel('subspace stength')
xlabel('subspace'); 
title('subspace strength colored motif (normalization')

%get the correlation between subspaces
[cstat,idx] = compareVectors(x,x,'corr',1);

figure; hold on; 
histogram(cstat)
ylabel('counts')
title('Pairwise similarity between subspace barcode across motifs');
xlabel('Rho (fisherZ)'); 

%project into full
xmat = idx2Matrix(cstat',idx');
figure; hold on; 
imagesc(xmat); colorbar; 
set(gca,'xlim',[0.5 size(xmat,1)],'ylim',[0.5 size(xmat,1)]); 
title('Pairwise similarity between subspace barcode across motifs');
xlabel('motif'); ylabel('motif'); 

%get the permutated distribution
rng('default')
n_perm = 1000;
cstat_perm = NaN(n_perm,numel(cstat));
for i = 1:n_perm
   perm_idx = arrayfun(@(n) randperm(size(x,1),size(x,1)),1:size(x,2),'UniformOutput',0);
   perm_idx = cat(1,perm_idx{:})';
   cstat_perm(i,:) = compareVectors(x,x(perm_idx),'sse',1);
end

%test the significance for each subspace, the difference between each
pval = NaN(1,size(cstat,1));
for i = 1:size(cstat,1)
    pval(i) = sum([cstat(i),cstat_perm(:,i)']<=cstat(i))/numel([cstat(i),cstat_perm(:,i)']);
end


end