function [score,idx] = bestLambda(cvl,maxflag,verbose)
%Camden - follows semendo et al., neuron 2019.
%best lambda of ridge is the largets value with mean value that lies within
%the 1 sem of the best model 
%Assumes cvl is ordered by fitting lambda from largest to smallest
%returns the cross-validated scores of that model, and it's indx

if nargin <2; maxflag=0; end %if 1, then returns the model with highest performance
if nargin <3; verbose=0; end 

% mean loss and standard error of the mean across folds.
cvLoss = [ mean(cvl); std(cvl)/sqrt(size(cvl,1)) ];

%get idx of bestLambda
y = 1-cvLoss(1,:);
e = cvLoss(2,:);

[~,best_idx] = max(y);
if maxflag==1
    idx = best_idx;
    score = 1-cvl(:,idx);
else
    idx = find(y>=(y(best_idx)-e(best_idx)),1,'first');
    score = 1-cvl(:,idx);
end

if verbose
   figure; hold on;
   errorbar(1:numel(y),y,e,'o--', 'Color','k')
   plot(idx,y(idx),'marker','x','color','r')
   xlabel('Lambda number')
   ylabel('Predictive performance')   
end %verbose

end %function end


