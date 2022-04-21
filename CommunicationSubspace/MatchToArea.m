function x_pad = MatchToArea(TruePairings,TrueArea,xPairings, xArea ,x)
%Camden MacDowell - timeless
%only works for xPairings<TruePairings
%x is cell

assert(numel(TruePairings)>=numel(xPairings),'True parings not greater than xPairings')

%if not a combination (depending on how processed)
if isvector(TruePairings)
    TruePairings = combvec(TruePairings,TruePairings)';
    TruePairings(TruePairings(:,2)-TruePairings(:,1)==0,:)=[];
end

if isvector(xPairings)
    xPairings = combvec(xPairings,xPairings)';
    xPairings(xPairings(:,2)-xPairings(:,1)==0,:)=[];
end
   
%match the xArea to TrueArea
x_pad = NaN(size(TruePairings,1),1);
for i = 1:size(TruePairings,1)
   %area 1 to area 2
   aTrue = find(ismember(TrueArea,TrueArea{TruePairings(i,1)}),1,'first');
   bTrue = find(ismember(TrueArea,TrueArea{TruePairings(i,2)}),1,'first');
   idx_true = intersect(find(TruePairings(:,1)==aTrue), find(TruePairings(:,2)==bTrue));
   
   aX = find(ismember(xArea,TrueArea{TruePairings(i,1)}),1,'first');
   bX = find(ismember(xArea,TrueArea{TruePairings(i,2)}),1,'first');
   if ~isempty(aX) && ~isempty(bX)
      idx = intersect(find(xPairings(:,1)==aX), find(xPairings(:,2)==bX));
   else
      idx = [];
   end
 
   %find the index within the pairings   
   if isempty(idx)    
       x_pad(idx_true) = NaN;
   else
       x_pad(idx_true) = x(idx);
   end
end

end %function end


