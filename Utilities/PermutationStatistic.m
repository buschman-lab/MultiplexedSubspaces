function perm_stat = PermutationStatistic(x,y,apply_func,n_shuf,usecell)

%permutes y and recomputes function f(x,y) 

if nargin < 5; usecell = 0; end
%permutation test
rng('default');

n = numel(y);
if usecell
    perm_stat = cell(1,n_shuf);
else
    perm_stat = NaN(1,n_shuf);
end
for i = 1:n_shuf
   temp_idx = randperm(n,n);
   if usecell
      perm_stat{i} = apply_func(x,y(temp_idx));
   else       
      perm_stat(i) = apply_func(x,y(temp_idx));
   end
end

% %OLD VERSION
% function perm_stat = PermutationStatistic(x,grp,apply_func,n_shuf,usecell)
% 
% %permutes y
% 
% if nargin < 5; usecell = 0; end
% %permutation test
% rng('default');
% 
% n = sum(grp==0);
% t = numel(grp);
% if usecell
%     perm_stat = cell(1,n_shuf);
% else
%     perm_stat = NaN(1,n_shuf);
% end
% for i = 1:n_shuf
%    temp_idx = randperm(t,n);
%    x1 = x(temp_idx);
%    x2 = x(~ismember([1:t],temp_idx));
%    if usecell
%       perm_stat{i} = apply_func(x1,x2);
%    else       
%       perm_stat(i) = apply_func(x1,x2);
%    end
% end


