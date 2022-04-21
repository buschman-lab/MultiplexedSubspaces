folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
%grab data
rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_mo\w*.mat'],0,{folder}); 

for cur_fn = 1:numel(fn)
    fprintf('\n\t working on file %d',cur_fn)
    load(fn{cur_fn});
    
    %parse activity per parent region 
    [area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general');

    %remove neurons that fire less than 0.5spike/sec on average across trials
    % [area_val, inactive_idx] = RemoveInactiveNeurons(area_val, 0.5/7.5);

    %clean up areas %third input is the min # of spikes to keep area
    [area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 
        
    
    for pairing = 1:size(paired_areas,1)
        fprintf('\n\t pairing %d',pairing)
        x = area_val{strcmp(area_label,area_label{paired_areas(pairing,1)}),:};
        y = area_val{strcmp(area_label,area_label{paired_areas(pairing,2)}),:};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');

        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);
     
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
        rng('default')
        dMaxShrink = .5:.01:1;
        lambda = GetRidgeLambda(dMaxShrink, x);

        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
            'LossMeasure', lossMeasure,'scale',false); 

        % Cross-validation routine.
        cvl_ridgetemp = crossval(cvFun, y, x, ...
              'KFold', cvNumFolds, ...
            'Options', cvOptions);

        %% Reduced Rank Regression XVal
        rng('default')
        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
            numDimsUsedForPrediction, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',false);

        % Cross-validation routine.
        cvl_rrrtemp = crossval(cvFun, y, x, ...
              'KFold', cvNumFolds, ...
            'Options', cvOptions);

        % Reduced Rank Regression Full with the optimal dimensionality
        d = ModelSelect([ mean(cvl_rrrtemp); std(cvl_rrrtemp)/sqrt(size(cvl_rrrtemp,1)) ], 1:size(cvl_rrrtemp,2));
        
        cvl_ridge{pairing} = cvl_ridgetemp;
        cvl_rrr{pairing} = cvl_rrrtemp;
        [~,rrr_B{pairing}] = ReducedRankRegress(y, x, d,'scale',false,'RIDGEINIT',false);
    end
    save(fn{cur_fn});
end