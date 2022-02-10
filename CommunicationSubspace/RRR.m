function [qOpt_source,qOpt_target,cvl_fa,rrr_B,cvl_rrr,cvl_ridge] = RRR(x,y)


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

%% Full model | Ridge regression
dMaxShrink = .5:.01:1;
lambda = GetRidgeLambda(dMaxShrink, x);

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
	'LossMeasure', lossMeasure);

% Cross-validation routine.
cvl_ridge = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

%% Reduced Rank Regression XVal
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
	numDimsUsedForPrediction, 'LossMeasure', lossMeasure);

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

% Reduced Rank Regression Full
[~, rrr_B] = ReducedRankRegress(y, x);


%% Factor analysis to determine dimensionalit of target neural population
q = 0:min([30,size(x,2)]);
cvLoss= CrossValFa(x, q, cvNumFolds, cvOptions);
qOpt_target = FactorAnalysisModelSelect(cvLoss, q);


%% Factor regression
% This finds the dominant dimensions in source activity and predicts the
% target with it
q = 0:min([size(x,2),30]);
qOpt_source = FactorAnalysisModelSelect( ...
	CrossValFa(x, q, cvNumFolds, cvOptions), ...
	q);

%here qOpt is the number of dimensions in the source

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@FactorRegress, Ytrain, Xtrain, Ytest, Xtest, ...
	numDimsUsedForPrediction, ...
	'LossMeasure', lossMeasure, 'qOpt', qOpt_source);

cvl_fa = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

% score = bestLambda(cvl_ridge,1); plot_rrrSummary(score,cvl_rrr,cvl_fa,0,1)
end %function end












