%see if the shared variance is just do to shared variance

load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1rrr_muaflag1_motif3.mat')

 for i = 7% 1:size(paired_areas,2)
    fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
%     idx = strcmp(area_label,area_label{paired_areas(i)});
%     x = cat(1,area_val{idx==0,:});
    x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
%     y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};
    y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');
    y = normalizeToBaseline(y,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);
    
    %get the trial averaged signal

%     %while here, go ahead and parse beta by area
%     grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
%     grp = cat(1,grp{:});
%     grouping{i} = grp;

    %subtract the psth
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);
    
    %add a neuron capturing the shared activity across neurons
    xx = cat(1,x,nanmean(x,1));
    yy = cat(1,y,nanmean(y,1));

    %concatentate across trials and pca
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
    yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';
   
    
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


    cvLoss = [1-mean(cvl_ridge); std(cvl_ridge)/sqrt(cvNumFolds) ];

    errorbar(1:numel(cvLoss(1,:)),cvLoss(1,:),cvLoss(2,:))

    [score,idx] = bestLambda(cvl_ridge);

    % Reduced Rank Regression XVal
    rng('default')
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));

    % Cross-validation routine.
    cvl_rrr = crossval(cvFun, y, x, ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);
    
    cvLoss = [1-mean(cvl_rrr); std(cvl_rrr)/sqrt(cvNumFolds) ];

    figure; hold on;
    errorbar(1:numel(cvLoss(1,:)),cvLoss(1,:),cvLoss(2,:))
    
    rng('default')
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));
    
    %% with trial mean added       
    rng('default')
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));
    
    cvl_rrr = crossval(cvFun, y, xx, ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);
    
    ylabel('performance'); xlabel('dimensions')
    cvLoss = [1-mean(cvl_rrr); std(cvl_rrr)/sqrt(cvNumFolds) ];
    
    errorbar(1:numel(cvLoss(1,:)),cvLoss(1,:),cvLoss(2,:))      
    
    % Reduced Rank Regression Full with the optimal dimensionality to get beta weights
    d = ModelSelect([ mean(cvl_rrr); std(cvl_rrr)/sqrt(size(cvl_rrr,1)) ], 1:size(cvl_rrr,2));
    [rrr_b,rrr_B,rrr_V] = ReducedRankRegress(y, xx, d,'scale',false,'RIDGEINIT',lambda(idx));   
    figure; t=tiledlayout(5,1);
    for jj = 1:5
        nexttile
        bar(rrr_B(:,jj));
        title(sprintf('dimension %d',jj))
    end
    

    
    [Q,R] = qr( rrr_B(:,1:5) );
    figure; t=tiledlayout(5,1);
    for jj = 1:5
        nexttile
        bar(Q(:,jj));
        title(sprintf('dimension %d',jj))
    end     
    dot(Q(:,1),Q(:,2))
        
    
    
        
        
end %subspace identification loop













% 
% 
% %% with trial mean
%     rng('default')
%     for jj = 1:size(xx,2)-1
%         xx(:,jj) = xx(randperm(size(xx,1),size(xx,1)),jj);        
%     end
%     for jj = 1:size(yy,2)-1
%         yy(:,jj) = yy(randperm(size(yy,1),size(yy,1)),jj);        
%     end
%    
%     rng('default')
%     dMaxShrink = .5:.01:1;
%     lambda = GetRidgeLambda(dMaxShrink, xx,'scale',false);
% 
%     cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
%         (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
%         'LossMeasure', lossMeasure,'scale',false); 
% 
%     % Cross-validation routine.
%     cvl_ridge = crossval(cvFun, y, xx, ...
%           'KFold', cvNumFolds, ...
%         'Options', cvOptions);
% 
% 
%     cvLoss = [1-mean(cvl_ridge); std(cvl_ridge)/sqrt(cvNumFolds) ];
% 
% %     hold on; errorbar(1:numel(cvLoss(1,:)),cvLoss(1,:),cvLoss(2,:))
% 
%     [score,idx] = bestLambda(cvl_ridge);    
%     