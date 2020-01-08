function new = BlockMatrix(orig,sz)

%splits matrix include equal sized, non overlapping blocks. 
%if given a blocked, matrix, will reconstuct into the original image. 
if iscell(orig)
    %to reconstruct the original image
    temp = cell(size(orig,1),1);
    for i = 1:size(orig,1)
       temp{i} = cat(2,orig{i,:});
    end
    new = cat(1,temp{:});
else
    [n,m]=size(orig);
    aa=1:sz:n;
    bb=1:sz:m;
    [ii,jj]=ndgrid(aa,bb);
    new=arrayfun(@(x,y) orig(x:x+sz-1,y:y+sz-1),ii,jj,'un',0);
end



