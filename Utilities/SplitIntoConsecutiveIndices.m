function y = SplitIntoConsecutiveIndices(x)

unique_val = unique(x);
%break each into a set of consecutive indices
indices = arrayfun(@(z) find(x==z),unique_val,'UniformOutput',0); 
consecutive_indices = cellfun(@(z) find(diff([false;[1;diff(z)]==1;false])~=0), indices,'UniformOutput',0);
consecutive_indices = cellfun(@(z) reshape(z', 2,[])', consecutive_indices,'UniformOutput',0);

%break into cell array of cells with the different consecutive indices
y= {};
for j = 1:numel(indices)
   x = indices{j};
   z = consecutive_indices{j};
   for i = 1:size(z,1)
       if i == 1
           y{j,i} = x(1:z(i,2)-1);
       elseif i == size(z,1)
           y{j,i} = x(z(i-1,2):end);
       else
           y{j,i} = x(z(i-1,2):z(i,2)-1);
       end
   end
end