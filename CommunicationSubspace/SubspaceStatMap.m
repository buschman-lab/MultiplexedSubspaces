function y = SubspaceStatMap(x,paired_areas,remove_diag)
%Camden MacDowell - timesless
%x is the paired_areaXnummotifs
if nargin <3; remove_diag=1; end

n = max(unique(paired_areas));
y = NaN(n,n,size(x,2));
for cur_m = 1:size(x,2)
   for i = 1:size(paired_areas,1)
       y(paired_areas(i,1),paired_areas(i,2),cur_m)= x(i,cur_m);
   end
end

if remove_diag
    for i = 1:size(y,3)
        temp = y(:,:,i);
        temp(1:1+size(y,1):end)=NaN;
        y(:,:,i) = temp;
    end
end

end %function end