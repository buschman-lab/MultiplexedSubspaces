function perm_stat = PermutationStatistic(x,grp,apply_func,n_shuf)
%permutation test
rng('default');

n = sum(grp==0);
t = numel(grp);
perm_stat = NaN(1,n_shuf);
for i = 1:n_shuf
   temp_idx = randperm(t,n);
   x1 = x(temp_idx);
   x2 = x(~ismember([1:t],temp_idx));
   perm_stat(i) = apply_func(x1,x2);
end

