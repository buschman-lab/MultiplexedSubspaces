function cvl_ridge = RidgePerm(x,y,permseed)

%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

if permseed >1
    rng(permseed)
    idx = randperm(size(x,3),size(x,3));
    x = x(:,:,idx);
end

%concatentate across trials 
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.

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


end %function end












