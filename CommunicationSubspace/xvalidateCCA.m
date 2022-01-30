function [r_test, deltaA, deltaB] = xvalidateCCA(x,y,n_folds)
%Camden MacDowell - timeless. Cross validates the Rsq for a given cca fit
%to test the generalizability. Also returns the variability in weights of
%the first CV (detlaA and deltaB)
if nargin <3; n_folds = 10; end

%get weights of full model
[a_full,b_full] = canoncorr(x,y);

rng('default'); 
c = cvpartition(size(x,1),'KFold',n_folds);
r_test = NaN(n_folds,min(size(x,2),size(y,2)));
deltaA = NaN(size(a_full,1),n_folds);
deltaB = NaN(size(b_full,1),n_folds);
for i = 1:n_folds
    %train
    [a,b] = canoncorr(x(c.training(i),:),y(c.training(i),:));
    %project testing data
    testU = x(c.test(i),:)*a;
    testV = y(c.test(i),:)*b;
    %get the r
    r_test(i,:) = arrayfun(@(n) corr(testU(:,n),testV(:,n)), 1:size(testU,2));
    deltaA(:,i)  = a_full(:,1)-a(:,1);
    deltaB(:,i)  = b_full(:,1)-b(:,1);
end
%convert to Rsq and avg
r_test = nanmean(r_test.^2);



end