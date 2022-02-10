function subspaceAsymmetry(best_idx,a,b,t)
%Camden MacDowell - timeless
%compares similarity between CVs in upper and lower triangle to see if
%multiple subspaces
  
%top triangle = first number >> second
idx = t{1};
pval=NaN(2,numel(best_idx));
for i = 1:numel(best_idx)
   within = []; between = [];
   temp = idx(:,best_idx{i});
   
   %within
   within(:,1) = abs(SubspaceConsistency(a{i}(:,diff(temp)<0),a{i}(:,diff(temp)<0),'corr',1));
   within(:,2) = abs(SubspaceConsistency(b{i}(:,diff(temp)<0),b{i}(:,diff(temp)<0),'corr',1));
   
   %between
   between(:,1) = abs(SubspaceConsistency(a{i}(:,diff(temp)<0),a{i}(:,diff(temp)>0),'corr',0)); 
   between(:,2) = abs(SubspaceConsistency(b{i}(:,diff(temp)<0),b{i}(:,diff(temp)>0),'corr',0)); 
   
   if sum(isnan(within(:)))==0 &&  sum(isnan(between(:)))==0
      pval(1,i) = ranksum(within(:,1),between(:,1),'tail','right');
      pval(2,i) = ranksum(within(:,2),between(:,2),'tail','right');       
   end
end


end %function 