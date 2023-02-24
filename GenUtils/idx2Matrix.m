function y = idx2Matrix(x,idx)
   n = max(unique(idx));
   y = NaN(n,n);
   for i = 1:size(idx,1)
       y(idx(i,1),idx(i,2))= x(i);
   end
end