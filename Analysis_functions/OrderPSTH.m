function x = OrderPSTH(x)
%Camden 
%x is a neuron x time spike raster. Organizes neurons by their peak FR

[~,idx] = max(x,[],2);
[~,idx] = sort(idx);
x= x(idx,:);

end