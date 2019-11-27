function d = dprime(x,y)

%x and y need to be observations x variables

d = (nanmean(x)+nanmean(y))/sqrt(0.5 * (nanvar(x)+nanvar(y)));

end