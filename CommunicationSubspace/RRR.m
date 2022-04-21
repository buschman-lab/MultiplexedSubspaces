function [qOpt_source,qOpt_target,cvl_fa,rrr_B,cvl_rrr,cvl_ridge,cvl_fa_target,rrr_Betaxvall,rrr_V] = RRR(x,y)


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
numDimsUsedForPrediction = 1:10;

% %pairwise correlation between neurons
% rho_across = corr(x,y);
% 
% rho_within = corr(x,x);
% rho_within(1:1+size(y,1):end)=NaN;

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
	d, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',1);

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

rrr_Betaxvall=[];

% Reduced Rank Regression Full with the optimal dimensionality
d = ModelSelect([ mean(cvl_rrr); std(cvl_rrr)/sqrt(size(cvl_rrr,1)) ], 1:size(cvl_rrr,2));
[~,rrr_B,rrr_V] = ReducedRankRegress(y, x, d,'scale',false,'RIDGEINIT',lambda(idx));

%% Factor analysis to determine dimensionality of target neural population
% rng('default')
% q = 0:min([30,size(y,2)]);
% cvLoss= CrossValFa(y, q, cvNumFolds, cvOptions);
% 
% %make sure enough dimensions were used to fully capture the shared variance
% %(i.e. at lease 99.5%)... if not, then use larger space (slower). 
% if sum(cvLoss<0.005)==1 %i.e. if only reaches it in the final iteration (which must ==1)
%     q = 0:min([50,size(y,2)]);
%     cvLoss= CrossValFa(y, q, cvNumFolds, cvOptions);
% end
% %save off the optimal number of dimensions
% qOpt_target = FactorAnalysisModelSelect(cvLoss, q);
% 
% %save off the full space
% cvl_fa_target = cvLoss; 
% 
% 
% %% Factor regression
% % This finds the dominant dimensions in source activity and predicts the
% % target with it
% rng('default')
% q = 0:min([size(x,2),30]);
% qOpt_source = FactorAnalysisModelSelect( ...
% 	CrossValFa(x, q, cvNumFolds, cvOptions), ...
% 	q);
% 
% cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
% 	(@FactorRegress, Ytrain, Xtrain, Ytest, Xtest, ...
% 	numDimsUsedForPrediction, ...
% 	'LossMeasure', lossMeasure, 'qOpt', qOpt_source);
% 
% cvl_fa = crossval(cvFun, y, x, ...
% 	  'KFold', cvNumFolds, ...
% 	'Options', cvOptions);

qOpt_source = [];
cvl_fa = [];
q = [];
cvl_fa_target = []; 
qOpt_target=[];
% score = bestLambda(cvl_ridge,1); plot_rrrSummary(score,cvl_rrr,cvl_fa,0,1)
end %function end












