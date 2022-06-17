function LeaveOneOutSubspaceFit(fullfitpath)
%Camden MacDowell - timeless
%refits the model while leaving out each area ot see it's contribution to
%said model
refitlambda = 0; %flag to refit lambda (takes ~16hrs) versus not refiting (<1hr)
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
end

load(fullfitpath,'paired_areas','area_label','area_val','cvl_rrr','cvl_ridge');

[savedir,name] = fileparts(fullfitpath);

% savepath = fullfile(savedir,[name,'_LOLambdaRefit.mat']);
savepath = fullfile(savedir,[name,'_LO0.mat']);

delta_rsq = NaN(size(paired_areas,2),size(paired_areas,2),30);
delta_rsq_rel = NaN(size(paired_areas,2),size(paired_areas,2),30);
delta_full = NaN(size(paired_areas,2),size(paired_areas,2),1);
for i = 1:size(paired_areas,2)
    fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
    idx = strcmp(area_label,area_label{paired_areas(i)});
    x = cat(1,area_val{idx==0,:});
    y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');
    y = normalizeToBaseline(y,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);

    %while here, go ahead and parse beta by area
    grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
    grp = cat(1,grp{:});
    unique_grp = unique(grp);
    
    if refitlambda==1
    else
        %subtract the psth
        x = x-nanmean(x,3);
        y = y-nanmean(y,3);

        %concatentate across trials and pca
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
        y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';    

        %hold lambda steady with full model   
        rng('default'); dMaxShrink = .5:.01:1;
        lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);
        [~,lambda_idx] = bestLambda(cvl_ridge{i});
        lambda = lambda(lambda_idx);
    end
    
    %Leave-one-out 
    for j = 1:numel(unique_grp)
        if refitlambda==1
            xloo = x(grp~=unique_grp(j),:,:);
            [~,~,~,loo_rrr] = RRR_full(xloo,y); %to refit lambda (50x as slow with minimal difference)
        else        
            xloo = x(:,grp~=unique_grp(j));
            [~,~,~,loo_rrr] = RRR_simple(xloo,y,lambda,1); 
            
            cvNumFolds = 10;
            lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
            cvOptions = statset('crossval'); % Initialize default options for cross-validation.

            %do the full model
            cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
                (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
                'LossMeasure', lossMeasure,'scale',false); 

            % Cross-validation routine.
            cvl_ridge_loo = crossval(cvFun, y, x, ...
                'KFold', cvNumFolds, ...
                'Options', cvOptions);
        end
        %get the difference
        temp = nanmean((1-cvl_rrr{i})-(1-loo_rrr))./nanmean((1-cvl_rrr{i})); %make relative to the total explainable variance of full model
        delta_rsq_rel(i,unique_grp(j),1:numel(temp)) = temp;
        temp = nanmean((1-cvl_rrr{i})-(1-loo_rrr));
        delta_rsq(i,unique_grp(j),1:numel(temp)) = temp;
        temp = nanmean((1-cvl_ridge{i}(:,lambda_idx))-(1-cvl_ridge_loo))./nanmean((1-cvl_ridge{i}(:,lambda_idx)));
        delta_full(i,unique_grp(j)) = temp;
    end
end %subspace identification loop

save(savepath,'delta_rsq','delta_rsq_rel','cvl_rrr','paired_areas','area_label','area_val','cvl_rrr','cvl_ridge','delta_full');


end %function end










