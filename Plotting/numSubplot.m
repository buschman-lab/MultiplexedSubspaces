function [num_rows, num_col]=numSubplot(n,r)

if nargin <2; r = 2; end %ratio of rows to columns;     

num_rows = floor(r*sqrt(n/r));
num_col = ceil(n/num_rows);

end
