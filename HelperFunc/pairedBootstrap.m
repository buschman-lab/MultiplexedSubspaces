function [xboot,stats] = pairedBootstrap(x,bfunc)
%camden macdowell - timeless
%performs a paired bootstrap for each column of x

rng('default');
xboot = NaN(1000,size(x,2));
for i = 1:1000
    idx = datasample(1:size(x,1),size(x,1));
    xboot(i,:) = bfunc(x(idx,:));    
end

stats.ci = cat(1,prctile(xboot,2.5),prctile(xboot,97.5));
stats.absolute = cat(1,prctile(xboot,0),prctile(xboot,100));
stats.mean = nanmean(xboot);

p = NaN(1,size(x,2));
for i = 1:size(x,2)
    p(i) = sum(xboot(:,i)<=0)/numel(xboot(:,i));
end
stats.p = p; %versus zero


end