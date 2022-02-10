function [e_ratio, e, pval, withinSEM] = errBetween(score,cvl)
%camden - timeless
%error between the full model (score) and optimal rrr model 
%cvl is the cross validations on the rrr model.  

d = ModelSelect([ mean(cvl); std(cvl)/sqrt(size(cvl,1)) ], 1:size(cvl,2));

%get the error between (not paired, since random cv initialization)
pval = ranksum(score,1-cvl(:,d));

e = nanmean(score-(1-cvl(:,d)));

e_ratio = e/(nanmean(score));

%does the model get within the sem of the full model
withinSEM = nanmean(1-cvl(:,d))>nanmean(score)-(std(score)/sqrt(numel(score)));

end %function end