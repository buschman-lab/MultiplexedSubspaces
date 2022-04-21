function shufidx = PermuteWithinArea(B,data,cur_motif,cur_area)

grp = data(cur_motif).grouping{cur_area};
if numel(grp)~=size(B,1) %working on the single area
   grp = ones(size(B,1),1);
end
ugrp = unique(grp);

shufidx = NaN(size(B));
for j = 1:size(B,2)
    tempidx = 1:numel(grp);
    for i = 1:numel(ugrp)
        temp = tempidx(grp==ugrp(i));
        %maintain sign of neuron
        s = B(grp==ugrp(i),j)>0;
        a = temp(s==1);
        temp(s==1) = a(randperm(numel(a),numel(a)));
        a = temp(s==0);
        temp(s==0) = a(randperm(numel(a),numel(a)));
        tempidx(grp==ugrp(i)) = temp;
    end
    shufidx(:,j) = tempidx;
end

end