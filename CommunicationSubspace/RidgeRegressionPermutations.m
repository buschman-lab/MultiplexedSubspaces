function RidgeRegressionPermutations(block,cur_rec,cur_motif)
%Camden MacDowell - timeless
%runs trial-shuffled permutation of the ridge regression to get baseline
%performance. 

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\RidgePermutations';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/RidgePermutations';
end

rng(block);
fprintf('Working on rec %d motif %d',cur_rec,cur_motif);
data = LoadSubspaceData('in_grouped');
fprintf('\nDone loading data');

areas = data{cur_rec}(cur_motif).area_label;
rsq = NaN(1,numel(areas));
rrr_B = cell(1,numel(areas));
rrr_V = cell(1,numel(areas));
for cur_a = 1:numel(areas)
    fprintf('\n Working on area %d of %d',cur_a,numel(areas));
    [x,y] = loadFunc(data,cur_rec,cur_motif,areas{cur_a}); 
    lambda = getLambda(data,cur_rec,cur_motif,areas{cur_a});
    x = conditionVar(x);     
    yperm = trialpermute(y);
    yperm = conditionVar(yperm);
    [rsq(cur_a),rrr_B{cur_a},rrr_V{cur_a}] = ridge_simple(x,yperm,lambda);
end
fn = [savedir filesep sprintf('shuf%drec%dm%d.mat',block,cur_rec,cur_motif)];
save(fn,'rsq','cur_rec','block','cur_motif','areas','rrr_B','rrr_V');
fprintf('\n Done saving data');
end %function 

function yperm = trialpermute(y)
    yperm = y(:,:,randperm(size(y,3),size(y,3)));
end


function lambda = getLambda(data,cur_rec,cur_motif,area_name) 
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;    
    area_idx = strcmp(area_label,area_name);
    
    %get the originally used lambda
    x = cat(1,area_val{area_idx==0});
    if size(x,3)>1
        x = x-nanmean(x,3);
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';        
    end
    dMaxShrink = .5:.01:1;
    lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);
    cvl_ridge = data{cur_rec}(cur_motif).cvl_ridge{area_idx};
    [~,idx] = bestLambda(cvl_ridge);
    lambda = lambda(idx);
    
end


function [x,y] = loadFunc(data,cur_rec,cur_motif,area_name)
    area_label = data{cur_rec}(cur_motif).area_label;
    area_val = data{cur_rec}(cur_motif).area_val;
    x = cat(1,area_val{strcmp(area_label,area_name)==0});
    y = area_val{strcmp(area_label,area_name)};

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');
    y = normalizeToBaseline(y,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);
    
    %subtract the psth
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);
end

function x = conditionVar(x)
    %concatentate across trials and pca
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
end

function [rsq,rrr_B,rrr_V] = ridge_simple(x,y,lambda)
% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:30;

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
    'LossMeasure', lossMeasure,'scale',false); 

% Cross-validation routine.
cvl_ridge = crossval(cvFun, y, x, ...
      'KFold', cvNumFolds, ...
    'Options', cvOptions);

rsq = 1-nanmean(cvl_ridge);

rng('default')
d = 1:min(numDimsUsedForPrediction(end),size(y,2));
% Reduced Rank Regression Full with the optimal dimensionality
[~,rrr_B,rrr_V] = ReducedRankRegress(y, x, d,'scale',false,'RIDGEINIT',lambda);


end %function end





