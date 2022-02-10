function [r_test, dThetaA, dThetaB, deltaA, deltaB] = xvalidateCCA(x,y,n_folds)
%Camden MacDowell - timeless. Cross validates the Rsq for a given cca fit
%to test the generalizability. Also returns the angle (dThetaA/dThetaB) between the full
%model weights in each region and the model. 
if nargin <3; n_folds = 10; end
rng('default')
c = cvpartition(size(x,1),'KFold',n_folds);
%remove neu with zero var in any fold. They are ignored in CCA
bad_idx = cell(n_folds,1);
for i = 1:n_folds
   bad_idx{i,1} = find(nanvar(x(c.training(i),:),[],1)<eps);
   bad_idx{i,2} = find(nanvar(y(c.training(i),:),[],1)<eps);   
end
x(:,unique([bad_idx{:,1}]))=[];
y(:,unique([bad_idx{:,2}]))=[];

%get weights of full model
[a_full,b_full,r,~,~,stats] = canoncorr(x,y);
% sig_idx = stats.p<0.05;

%run cross validations
r_test = NaN(n_folds,min(size(x,2),size(y,2)));
dThetaA = NaN(n_folds,1);
dThetaB = NaN(n_folds,1);
deltaA = NaN(size(a_full,1),n_folds);
deltaB = NaN(size(b_full,1),n_folds);
for i = 1:n_folds   
    %train
    [a,b] = canoncorr(x(c.training(i),:),y(c.training(i),:));
    %project testing data
    testU = x(c.test(i),:)*a;
    testV = y(c.test(i),:)*b;  
    
    %get the r
    r_test(i,1:size(testU,2)) = arrayfun(@(n) corr(testU(:,n),testV(:,n)), 1:size(testU,2));
    dThetaA(i) = AngleBetweenWeights(a_full(:,1),a(:,1),'none'); %just get the first variate, as this will have the largest magnitude
    dThetaB(i) = AngleBetweenWeights(b_full(:,1),b(:,1),'none');
    deltaA(:,i)  = a_full(:,1)-a(:,1); %just get the first variate, as this will have the largest magnitude
    deltaB(:,i)  = b_full(:,1)-b(:,1);
end
%convert to Rsq and avg
r_test = nanmean(r_test.^2);



end