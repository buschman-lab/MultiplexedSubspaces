function [rrr_b,rrr_B,rrr_V,cvl_rrr,cvl_ridge,cvl_rrr_unique,contribution,rrr_Betaxvall] = RRR_full(x,y,grp)
if nargin <3; grp=[]; end

%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%concatentate across trials and pca
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:30;

%% Full model | Ridge regression
rng('default')
dMaxShrink = .5:.01:1;
lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
	'LossMeasure', lossMeasure,'scale',false); 

% Cross-validation routine.
cvl_ridge = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

[~,idx] = bestLambda(cvl_ridge);

%% Reduced Rank Regression XVal
rng('default')
d = 1:min(numDimsUsedForPrediction(end),size(y,2));
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
	d, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

rrr_Betaxvall=[];

% Reduced Rank Regression Full with the optimal dimensionality to get beta weights
d = ModelSelect([ mean(cvl_rrr); std(cvl_rrr)/sqrt(size(cvl_rrr,1)) ], 1:size(cvl_rrr,2));
[rrr_b,rrr_B,rrr_V] = ReducedRankRegress(y, x, d,'scale',false,'RIDGEINIT',lambda(idx));

%% get the unique contribution of each area
% if ~isempty(grp)
%     rng('default')
%     unique_grp = unique(grp);       
%     cvl_rrr_unique = NaN(cvNumFolds,numel(unique_grp));
%     for i = 1:numel(unique_grp)
%         temp = x(:,grp==unique_grp(i));
%         %shuffle all timepoints in a unique manner
%         temp = arrayfun(@(n) temp(randperm(size(temp,1),size(temp,1)),n), 1:size(temp,2),'UniformOutput',0);
%         temp = cat(2,temp{:});
%         xx=x;
%         xx(:,grp==unique_grp(i))=temp;
%         cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
%             (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
%             d, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));        
%         % Cross-validation routine.
%         cvl_rrr_unique(:,i) = crossval(cvFun, y, xx, ...
%               'KFold', cvNumFolds, ...
%             'Options', cvOptions);        
%         
%     end    
% end
% 
% full_rsq = (1-nanmean(cvl_rrr(:,d)));
% unique_rsq = (1-nanmean(cvl_rrr_unique));
% contribution = (full_rsq-unique_rsq)/full_rsq;
contribution=[];
unique_rsq=[];
full_rsq=[];
cvl_rrr_unique=[];
end %function end












