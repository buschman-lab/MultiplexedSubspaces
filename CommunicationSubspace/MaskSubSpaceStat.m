function x = MaskSubSpaceStat(x,idx,id)
%camden 
%x is a paired_idx x motif matrix. idx is a 1x motif cell array
%id (def==0) is the number to consider 'bad' in the index (too allow
%for masking by good/bad or diff strengths). 

if nargin <3; id =0 ; end

for i = 1:numel(idx)
    x(idx{i}==id,i)=NaN;
end


end %function end