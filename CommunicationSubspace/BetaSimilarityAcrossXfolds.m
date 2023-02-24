function r = BetaSimilarityAcrossXfolds(cur_rec)
%Camden MacDowell - timeless
%Code to analyze the results is included at the bottom. 
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

data = LoadSubspaceData('in_grouped');
data = data{cur_rec};
%Camden - timeless
ndim = 10;
nxval = 10;
rng('default')
area_label = data(1).area_label;
%loop through each area
r = NaN(numel(area_label),size(data,2),ndim);
for cur_area = 1:numel(area_label) 
    fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));
    tic
    for cur_motif = 1:size(data,2)
        x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});
        y = data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==1};

        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');
        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);
        %remove the psth
        x = x-nanmean(x,3);
        y = y-nanmean(y,3);
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';             
        y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';        
        
        %10 fold cross validation        
        cvp = cvpartition(size(x,1),'kfold',nxval);
        
        %get the lambda (using full rec as done with initial fitting)
        dMaxShrink = .5:.01:1;
        lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);        
        [~,idx] = bestLambda(data(cur_motif).cvl_ridge{cur_area});
        
        % Cross validated Betas
        B = cell(1,nxval); V = cell(1,nxval);
        for i = 1:nxval
            [~,B{i},V{i}] = ReducedRankRegress(y(cvp.training(i),:), x(cvp.training(i),:), ndim,'scale',false,'RIDGEINIT',lambda(idx));               
        end
        
        %get correlation by dimension (note that this will do ALL dim)
        for i = 1:size(B{1},2)
            x = cellfun(@(x) x(:,i),B,'UniformOutput',0); 
            x = cat(2,x{:});
            rho = corr(x);
            rho(triu(ones(nxval),0)==1) = NaN;
            r(cur_area,cur_motif,i) = nanmean(rho(:));            
        end        
    end
    toc
end %area

if ispc
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CrossValidationSimilarity\';
else
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CrossValidationSimilarity/';
end
save([savedir, sprintf('rec%d',cur_rec)],'r');

end %function 


%% code to analysze the results
fp = fig_params_cortdynamics;
fn = GrabFiles('\w*.mat',0,{pwd});
data = cellfun(@(x) load(x), fn);
x = {};
ndim = 50; 
for i = 1:numel(data)
    temp = data(i).r; 
    temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));    
    x{i} = temp(:,1:ndim);
end
x = cat(1,x{:});
temp = fisherZ(x); 

figure; hold on; 
[a,stats] = pairedBootstrap(temp,@nanmean);
stats.ci = fisherInverse(stats.ci);
stats.mean = fisherInverse(stats.mean);
errorbar(1:50,stats.mean,stats.mean-stats.ci(1,:),stats.ci(2,:)-stats.mean,'marker','o',...
    'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','none','linewidth',1.5,'color',[0.5 0.5 0.5])
% CompareViolins(a',fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',repmat({[0.15 0.15 0.15]},1,size(x,2)),'distWidth',0.75,'connectline',[0.75 0.75 0.75]);
xlim([0 size(x,2)+1])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 10 4],[2 10 20 10])
set(gca,'xtick',0:5:50,'xticklabel',0:5:50);
xlabel('Subspace dimension');
ylabel('Rho')
ylim([0 1])
title({'Similarity of Fitted B','weights across dimensions'},'fontweight','normal');
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'XvalidationRho','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SimilarityXValidationBetas',0); close all

%% stats for text
a = temp(:,1:10);
[~,stats] = pairedBootstrap(a(:),@nanmean);
fisherInverse(stats.mean)
fisherInverse(stats.ci)















