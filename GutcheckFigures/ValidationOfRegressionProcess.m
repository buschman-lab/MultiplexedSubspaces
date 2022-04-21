function ValidationOfRegressionProcess(pairing)
%camden macdowell timeless

load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1RRR_muaflag1_motif3.mat',...
    'paired_areas','st_norm','area_val','area_label','motif_onset','motif')

% pairing=60;
fprintf('\n\tWorking on subspace pairing %d of %d',pairing,size(paired_areas,1));
x = area_val{strcmp(area_label,area_label{paired_areas(pairing,1)}),:};
% idx = strcmp(area_label,area_label{paired_areas(pairing,2)});
% x = cat(1,area_val{idx==0,:});
y = area_val{strcmp(area_label,area_label{paired_areas(pairing,2)}),:};

%normalize to baseline
x = normalizeToBaseline(x,[1:2],'mean');
y = normalizeToBaseline(y,[1:2],'mean');

%use post stimulus
x = x(:,3:end,:);
y = y(:,3:end,:);

%subtract the psth
x_avg = nanmean(x,3);
y_avg = nanmean(y,3);
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%concatentate across trials and pca
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

% Show X and Y
figure('position',[469 67 1073 894]); 
subplot(1,2,1); 
imagesc(x,[prctile(x(:),1),prctile(x(:),95)]); c=colorbar;
ylabel(c, {'Firing rate','(normalized and mean-substracted'})
title([area_label{paired_areas(pairing,1)} ' residuals | motif', num2str(motif)],'FontWeight','normal')
ylabel('trials x timepoints'); xlabel('neurons/MUAs');

subplot(1,2,2);
imagesc(y,[prctile(y(:),1),prctile(y(:),95)]); c=colorbar;
ylabel(c, {'Firing rate','(normalized and mean-substracted'})
title([area_label{paired_areas(pairing,2)} ' residuals | motif', num2str(motif)],'FontWeight','normal')
ylabel('trials x timepoints'); xlabel('neurons/MUAs'); colormap magma
sgtitle('Trial-to-trial variability in neural activity');

% Show X and Y
figure; hold on; 
subplot(1,2,1); 
imagesc(x_avg,[prctile(x_avg(:),1),prctile(x_avg(:),99)]); c=colorbar;
ylabel(c, {'Firing rate','(baseline normalized)'})
title([area_label{paired_areas(pairing,1)} ' Trial-average | motif', num2str(motif)],'FontWeight','normal')
xlabel('time from onset (ms)'); ylabel('neurons/MUAs');
set(gca,'xtick',[1,5,10],'XTickLabel',[0,4,9]*134);

subplot(1,2,2);
imagesc(y_avg,[prctile(y_avg(:),1),prctile(y_avg(:),99)]); c=colorbar;
ylabel(c, {'Firing rate','(baseline normalized)'})
title([area_label{paired_areas(pairing,2)} ' Trial-average | motif', num2str(motif)],'FontWeight','normal')
xlabel('time from onset (ms)'); ylabel('neurons/MUAs'); colormap magma
sgtitle('Trial-average response');
set(gca,'xtick',[1,5,10],'XTickLabel',[0,4,9]*134);


%% Get ridge regression lambda
% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:10;
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
lambda_opts=lambda(idx);

%fit the full model (where no penality will be the best fit)
B = RidgeRegress(y,x,lambda);
[loss_ridge, ~] = RegressPredict(y,x, B,'lossMeasure',lossMeasure);

%plot with the cross validation: basically asks how well the fit
%generalizes with increasing regularization
cvLoss = [ mean(cvl_ridge); std(cvl_ridge)/sqrt(cvNumFolds) ];
figure; hold on; 
errorbar(lambda, 1-cvLoss(1,:),cvLoss(2,:),'o--','color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'Markersize',10)
errorbar(lambda(idx), 1-cvLoss(1,idx),cvLoss(2,idx),'o--','color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'Markersize',10)
set(gca,'xscale','linear')
%plot the performance per lambda of fully trained model
plot(lambda,1-loss_ridge,'color','k','linewidth',2)
plot(lambda(idx),1-loss_ridge(idx),'marker','o','color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'Markersize',10)

title('Identifying Ridge \lambda','fontweight','normal');
xlabel('\lambda')
ylabel('Performance');
legend('Cross-validated','Optimal \lambda','fully-trained','Opt \lambda on fully trained')

% figure;
% Bridge = RidgeRegress(y,x,lambda(idx));
% bar(Bridge)
%% full model and full Ridge Regression should eventually meet
B = RidgeRegress(y,x,lambda_opts);
[loss_ridge_best, ~] = RegressPredict(y,x, B,'lossMeasure',lossMeasure);
B = RidgeRegress(y,x,0);
[loss_ridge_none, ~] = RegressPredict(y,x, B,'lossMeasure',lossMeasure);

n = min(size(x,2),size(y,2));
loss_rrr = NaN(1,n);
for d = 1:n
    B = ReducedRankRegress(y, x, d);
    [loss_rrr(d), ~] = RegressPredict(y, x, B,'lossMeasure','NSE');
end

figure; hold on; 
plot(1:numel(loss_rrr),1-loss_rrr,'marker','none','linestyle','-','color',[0.8 0.1 0.1],'linewidth',2)
plot([0,numel(loss_rrr)],[1-loss_ridge_best,1-loss_ridge_best],'linestyle',':','color',[0.1 0.1 0.8],'linewidth',1.5);
plot([0,numel(loss_rrr)],[1-loss_ridge_none,1-loss_ridge_none],'linestyle','--','color',[0.1 0.1 0.8],'linewidth',1.5);
title('Full Rank RRR Matches Ridge Performance (No xvalidation)','fontweight','normal')
ylabel('Performance');
xlabel('Rank (i.e., dimensionality)')
legend('RRR','Ridge opt \lambda','Ridge no \lambda','Location','southeast')


%% Reduced Rank Regression XVal same procedure
% RRR (no regularization)
n = min(size(x,2),size(y,2));
rng('default')
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
    1:n, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',0);

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
      'KFold', 10, ...
    'Options', cvOptions);
cvLoss = [ mean(cvl_rrr); std(cvl_rrr)/sqrt(cvNumFolds) ];
cvLoss_ridge = [ mean(cvl_ridge); std(cvl_ridge)/sqrt(cvNumFolds) ];

%repeat with regularization
rng('default')
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
    1:n, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda_opts);

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
      'KFold', 10, ...
    'Options', cvOptions);
cvLoss_rrr_ridge = [ mean(cvl_rrr); std(cvl_rrr)/sqrt(cvNumFolds) ];

% [B,B_,V] = ReducedRankRegress(y,x,1,'scale',false,'RIDGEINIT',lambda_opts);
% [rrr_best_loss, ~] = RegressPredict(y,x, B,'lossMeasure',lossMeasure);

figure('position',[681         559        1024         420]); hold on; 
errorbar(1:n, 1-cvLoss(1,:),cvLoss(2,:),'o--','color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'Markersize',5)
errorbar(1:n, 1-cvLoss_rrr_ridge(1,:),cvLoss_rrr_ridge(2,:),'o--','color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',5)
plot([1,n], [1-min(nanmean(cvl_ridge)),1-min(nanmean(cvl_ridge))],'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'Markersize',10,'linestyle','-','linewidth',1.5)
% plot([1,n], [1-cvLoss_ridge(1,idx),1-cvLoss_ridge(1,idx)],'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'Markersize',10,'linestyle',':','linewidth',1.5)
% plot([1,n], [1-cvLoss_ridge(1,end),1-cvLoss_ridge(1,end)],'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'Markersize',10,'linestyle','--','linewidth',1.5)
ylabel('Performance');
xlabel('Rank (i.e., dimensionality)')
legend('RRR','regularized RRR','Ridge peak \lambda','Ridge opt \lambda','Ridge no \lambda','Location','east')
title('vanilla RRR overfits resulting in poor xvalidation | regularized RRR matches expected performance','fontweight','normal')

%
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspacesTemp';
saveCurFigs(get(groot, 'Children'),{'-dpng'},['ValidationOrRegressionExamples',sprintf('pairing %d motif %d','', pairing, motif)],savedir,0); close all

end
























