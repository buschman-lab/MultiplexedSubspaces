function [a_full,b_full,U_full,V_full,rsq_xv,pvals,t,best_idx,aTheta_xv,bTheta_xv] = significantCVs_perTimePoint(x,y,alpha,verbose)
%Camden MacDowell - timeless
%follows methods described in: https://www.nature.com/articles/s41467-020-17902-1#Sec13
%"single-trial cross-area neural population dynamics during long-term skill
%learning. Veuthey et al., Nat Comm. 2020. 
%performs cca and determine generalizability with cross validation and
%significant cannical variables with permutation
%inputs: 
%x and y are two tensors neurons x timepoits x trials. 
%alpha is the significance lvl

if nargin <3; alpha = 0.05; end 
if nargin <4; verbose = 1; end

%get all pairs of timepoints between the two regions
t = combvec(1:size(x,2),1:size(y,2));

rsq_xv = cell(1,size(t,2));
deltaA = cell(1,size(t,2));
deltaB = cell(1,size(t,2));
for i = 1:size(t,2)
   xx = squeeze(x(:,t(1,i),:))';
   yy = squeeze(y(:,t(2,i),:))';
   [rsq_xv{i},deltaA{i}, deltaB{i}] = xvalidateCCA(xx,yy,10);
end
%adjust in case unequal lengths
n = min(cellfun(@numel,rsq_xv)); 
rsq_xv = cellfun(@(x) x(1:n), rsq_xv,'UniformOutput',0);
rsq_xv = cat(1,rsq_xv{:});

%figure example | plot the xvalidated xcorrellelogram
if verbose
    figure; hold on; 
    rsq_mat = NaN(size(x,2),size(y,2));
    for i = 1:size(t,2)
       rsq_mat(t(1,i),t(2,i)) = rsq_xv(i,1); 
    end
    imagesc(rsq_mat); colorbar; colormap magma
    xlim([0.5 size(rsq_mat,1)]); ylim([0.5 size(rsq_mat,1)]);
    set(gca,'xtick',[1:5:size(rsq_mat,1)],'xticklabel',[-5:5:15]);
    set(gca,'ytick',[1:5:size(rsq_mat,1)],'yticklabel',[-5:5:15]);
    xlabel('area 1'); ylabel('area 2')
end

%% Significance and stability testing versus trial-permuted data
%rank first CVs by their timepoints and use strongest for permutation test
[r_sort,idx] = sort(rsq_xv(:,1),'descend');

%trial-wise permutation of strongest CV
xx_best = squeeze(x(:,t(1,idx(1)),:))';
yy_best = squeeze(y(:,t(2,idx(1)),:))';
rng('default')
n_perm = 1000;

%preallocate permutations since xvalidate will reset rng
yy_perm = arrayfun(@(n) yy_best(randperm(size(yy_best,1),size(yy_best,1)),:),1:n_perm,'UniformOutput',0);
[rsq_xv_perm,dap,dbp] = arrayfun(@(n) xvalidateCCA(xx_best,yy_perm{n},10),1:n_perm,'UniformOutput',0);
rsq_xv_perm = arrayfun(@(n) rsq_xv_perm{n}(1),1:n_perm);

%test timepoint significance and stability
pval_r = arrayfun(@(x) sum([rsq_xv_perm,x]>=x)/numel([rsq_xv_perm,x]),r_sort);
 
a = arrayfun(@(n) nanmean(deltaA{n},'all'), idx);
a_perm = cellfun(@(x) nanmean(x,'all'),dap);
pval_a = arrayfun(@(x) sum([a_perm,x]<=x)/numel([a_perm,x]),a);

b = arrayfun(@(n) nanmean(deltaB{n},'all'), idx);
b_perm = cellfun(@(x) nanmean(x,'all'),dbp);
pval_b = arrayfun(@(x) sum([b_perm,x]<=x)/numel([b_perm,x]),b);

%1st CV for each time is independent so take any that meet all criteria 
sig_idx = find(sum(cat(2,pval_r,pval_a,pval_b)<alpha,2)==3);

%visualize the map of significant CV
if verbose    
    figure; hold on; 
    rsq_mat = NaN(size(x,2),size(y,2));
    for i = 1:numel(sig_idx)
       rsq_mat(t(1,idx(sig_idx(i))),t(2,idx(sig_idx(i)))) = r_sort(sig_idx(i),1); 
    end
    imagesc(rsq_mat); colorbar; colormap magma
    xlim([0.5 size(rsq_mat,1)]); ylim([0.5 size(rsq_mat,1)]);
    set(gca,'xtick',[1:5:size(rsq_mat,1)],'xticklabel',[-5:5:15]);
    set(gca,'ytick',[1:5:size(rsq_mat,1)],'yticklabel',[-5:5:15]);
    xlabel('area 1'); ylabel('area 2')
end

%run CCA on the full model at each significant timepoint
a_full = NaN(size(x,1),numel(sig_idx));
b_full = NaN(size(y,1),numel(sig_idx));
U_full = NaN(size(x,3),numel(sig_idx));
V_full = NaN(size(x,3),numel(sig_idx));
for i = 1:numel(sig_idx)
   xx = squeeze(x(:,t(1,idx(sig_idx(i))),:))';
   yy = squeeze(y(:,t(2,idx(sig_idx(i))),:))';    
   [aa,bb,~,uu,vv] = canoncorr(xx,yy);
   a_full(:,i) = aa(:,1);
   b_full(:,i) = bb(:,1);
   U_full(:,i) = uu(:,1);
   V_full(:,i) = vv(:,1);
end

%return significant cvs
rsq_xv = r_sort(sig_idx);
aTheta_xv = a(sig_idx);
bTheta_xv = b(sig_idx);
best_idx = idx(sig_idx);
pvals = cat(2,pval_r,pval_a,pval_b);
pvals = pvals(sig_idx,:);

%in future, you could add:
%if any timepoints are significant, then test the subsequent CVs at that
%timepoint

end %function end















