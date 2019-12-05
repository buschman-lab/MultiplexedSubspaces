function y = SplitIntoSubarray(y,len)

remainder = mod(size(y,1),len);
y(end-remainder+1:end,:) = [];
y = reshape(y, len, size(y,1)/len);

end