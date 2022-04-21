function TrialNumberAndNeuronNumberImpact()
%Camden - timeless
%Gut check to plot the impact of number of neurons on strength of subspace
%and number of motif occurances on strength of subspace

%number of neurons across motifs 
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA';
%grab data
rec_name = 'Mouse 331 Recording 1';
[fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_mo\w*.mat'],0,{folder}); 
fn(2)=[]; fn(end)=[];
data_rrr = cellfun(@(x) load(x),fn);
 
[~,~,~,n] = localDimensionality(data_rrr(1));
pairing = data_rrr(1).paired_areas;
%number of predictor neurons
n_pred = arrayfun(@(nn) n(nn), pairing(:,1));
%number of target neurons
n_targ = arrayfun(@(nn) n(nn), pairing(:,2));
score=[];
n_trial = [];
for cur_fit = 1:size(data_rrr,2)
    score(:,cur_fit) = cellfun(@(x)nanmean(bestLambda(x)), data_rrr(cur_fit).cvl_ridge,'UniformOutput',1);
    n_trial(cur_fit) = size(data_rrr(cur_fit).x,1);
end
n_pred(score(:,1)==1)=[];
n_targ(score(:,1)==1)=[];
score(score(:,1)==1,:)=[];

n_trial = repmat(n_trial,size(score,1),1);
n_pred = repmat(n_pred,size(score,2),1);
n_targ = repmat(n_targ,size(score,2),1);
score = score(:);
n_trial = n_trial(:);

%number of predictive neurons across motifs
figure; hold on; 
lm = fitlm(n_pred,score);
plot(lm);
[rho,p] = corr(n_pred,score,'type','Spearman');
title(sprintf('Across motifs spearman rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# predictive neurons')

%number of target neurons across motifs
figure; hold on; 
lm = fitlm(n_targ,score);
plot(lm);
[rho,p] = corr(n_targ,score,'type','Spearman');
title(sprintf('Across motifs spearman rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# target neurons')

%number of trials across motifs
figure; hold on; 
lm = fitlm(n_trial,score);
plot(lm);
[rho,p] = corr(n_trial,score,'type','Spearman');
title(sprintf('Across motifs spearman rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# trials')

%% not do subsampling of neurons within a motif
load(fn{3},'paired_areas','st_norm','area_val','area_label','motif_onset','motif')

n_perm=1;
paired_areas(paired_areas(:,1)-paired_areas(:,2)==0,:)=[];
targ_sweep_all = cell(1,size(paired_areas,1));
pred_sweep_all = cell(1,size(paired_areas,1));
targ_num = cell(1,size(paired_areas,1));
pred_num = cell(1,size(paired_areas,1));
for pairing = 1:size(paired_areas,1)
    fprintf('\n\tWorking on subspace pairing %d of %d',pairing,size(paired_areas,1));
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
    
    %range of # neurons
    pred = round(linspace(round(0.25*size(x,2)),size(x,2),4));
    targ = round(linspace(round(0.25*size(y,2)),size(y,2),4));
    
    %predicting    
    rng('default')
    pred_sweep = NaN(numel(pred),n_perm,2);
    for ii = 1:numel(pred)
        for jj = 1:n_perm %5 random subsamplings
           xx = x(:,randperm(size(x,2),pred(ii)));
           dMaxShrink = .5:.01:1;
           lambda = GetRidgeLambda(dMaxShrink, xx,'scale',false);

           cvOptions = statset('crossval');
           cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
               (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
                'LossMeasure', 'NSE','scale',false); 
           % Cross-validation routine.
           cvl_ridge = crossval(cvFun, y, xx, ...
                  'KFold', 10, ...
                'Options', cvOptions);


           [~,idx] = bestLambda(cvl_ridge);
           loss = [ 1-mean(cvl_ridge); std(cvl_ridge)/sqrt(10) ];
           pred_sweep(ii,jj,:) = loss(:,idx);
        end
    end
    
    %target neurons    
    rng('default')
    targ_sweep = NaN(numel(targ),n_perm,2);
    for ii = 1:numel(targ)
        for jj = 1:n_perm %5 random subsamplings
           yy = y(:,randperm(size(y,2),targ(ii)));
           dMaxShrink = .5:.01:1;
           lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);

           cvOptions = statset('crossval');
           cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
               (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
                'LossMeasure', 'NSE','scale',false); 
           % Cross-validation routine.
           cvl_ridge = crossval(cvFun, yy, x, ...
                  'KFold', 10, ...
                'Options', cvOptions);


           [~,idx] = bestLambda(cvl_ridge);
           loss = [ 1-mean(cvl_ridge); std(cvl_ridge)/sqrt(10) ];
           targ_sweep(ii,jj,:) = loss(:,idx);
        end
    end       
    
    %average across permutations
    targ_num{pairing} = targ;
    pred_num{pairing} = pred;
    targ_sweep_all{pairing} = squeeze(nanmean(targ_sweep,2));
    pred_sweep_all{pairing} = squeeze(nanmean(pred_sweep,2));   
end

%# of predictors
a = cat(1,pred_num{:});
b = cellfun(@(x) x(:,1),pred_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on; 
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Predictor # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# predictive neurons')

a = cat(2,pred_num{:})';
b = cellfun(@(x) x(:,2),pred_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on; 
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Predictor # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('instability (SEM of NSE)')
xlabel('# predictive neurons')

%# of targets
a = cat(2,targ_num{:})';
b = cellfun(@(x) x(:,1),targ_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on; 
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Target # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# Target neurons')

a = cat(2,targ_num{:})';
b = cellfun(@(x) x(:,2),targ_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on; 
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Target # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('instability (SEM of NSE)')
xlabel('# Target neurons')

%% impact number of trials

n_perm=1;
trial_num = cell(1,size(paired_areas,1));
trial_sweep_all = cell(1,size(paired_areas,1));
for pairing = 1:size(paired_areas,1)
    fprintf('\n\tWorking on subspace pairing %d of %d',pairing,size(paired_areas,1));
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

    %range of trial numbers
    trin = round(linspace(round(0.25*size(x,3)),size(x,3),4));
    
    rng('default');
    for ii = 1:numel(trin)
        for jj = 1:n_perm %5 random subsamplings            
            rand_idx = randperm(size(x,3),trin(ii));
            
            xx = x(:,:,rand_idx);
            yy = y(:,:,rand_idx);
            
            %concatentate across trials and pca
            xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';
            yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';

            dMaxShrink = .5:.01:1;
            lambda = GetRidgeLambda(dMaxShrink, xx,'scale',false);

            cvOptions = statset('crossval');
            cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
               (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
                'LossMeasure', 'NSE','scale',false); 
            % Cross-validation routine.
            cvl_ridge = crossval(cvFun, yy, xx, ...
                  'KFold', 10, ...
                'Options', cvOptions);


            [~,idx] = bestLambda(cvl_ridge);
            loss = [ 1-mean(cvl_ridge); std(cvl_ridge)/sqrt(10) ];
            trial_sweep(ii,jj,:) = loss(:,idx);            
        end
    end
    %average across permutations
    trial_num{pairing} = trin;
    trial_sweep_all{pairing} = squeeze(nanmean(trial_sweep,2));
end

%# of targets
a = cat(2,trial_num{:})';
b = cellfun(@(x) x(:,1),trial_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on;
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Trial # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('Performance (NSE)')
xlabel('# Trials')

a = cat(2,trial_num{:})';
b = cellfun(@(x) x(:,2),trial_sweep_all,'UniformOutput',0);
b = cat(2,b{:})';
a = a(:);
b = b(:);
figure; hold on; 
lm = fitlm(a,b);
plot(lm);
[rho,p] = corr(a,b,'type','Spearman');
title(sprintf('Impact of Trial # rho = %0.2f pval = %0.2g',rho,p),'fontweight','normal')
ylabel('instability (SEM of NSE)')
xlabel('# Trials')

save_dir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\TrialAndNeuronNumberGutcheck'; 
save([save_dir filesep 'TrialNeuNumberGutcheck.mat']);
saveCurFigs(get(groot, 'Children'),{'-dpng'},'TrialNeuNumberGutcheck',save_dir,0); %close all;












    