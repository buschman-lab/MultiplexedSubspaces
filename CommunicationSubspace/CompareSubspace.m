function CompareSubspace(B1,B2,x1,x2,y1,y2)
%Camden MacDowell - timeless
%

load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1RRR_muaflag1_motif3.mat',...
    'area_label','paired_areas','area_val')
i=2;
x = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
y = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};

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


load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA\Mouse 331 Recording 1RRR_muaflag1_motif4.mat',...
    'area_label','paired_areas','area_val')
i=2;
xx = area_val{strcmp(area_label,area_label{paired_areas(i,1)}),:};
yy = area_val{strcmp(area_label,area_label{paired_areas(i,2)}),:};

%normalize to baseline
xx = normalizeToBaseline(xx,[1:2],'mean');
yy = normalizeToBaseline(yy,[1:2],'mean');

%use post stimulus
xx = xx(:,3:end,:);
yy = yy(:,3:end,:);

%subtract the psth
xx = xx-nanmean(xx,3);
yy = yy-nanmean(yy,3);

%concatentate across trials and pca
xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';
yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';


%get the explained variance when predicting something else versus some
%measure of self (cross validations? Withheld, something like that)
motifs = [1,2];
%grab betas for motif 1 and motifs 2

%

rng('default')
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
	numDimsUsedForPrediction, 'LossMeasure', lossMeasure);

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);

end

