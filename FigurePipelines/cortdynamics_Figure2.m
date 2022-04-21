%Figure 2 
%Camden MacDowell


%% Load subspace data
data = LoadSubspaceData('in_grouped');
dataout = LoadSubspaceData('out_grouped');

%% Plot an example predictors and predicted

x = cat(1,data{1}(5).area_val{1:7});
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';



%%if you ever want to rerun the SharedBeta
% rsq = cellfun(@(x) SharedBetas_InternalValidation(x,0), data,'UniformOutput',0);
% rsqout = cellfun(@(x) SharedBetas_InternalValidation(x,1),dataout,'UniformOutput',0);
% adjust for the fact no RSP and bfd in rec 3 and 4
% temp = NaN(8,14,10);
% temp([1:3,5,7,8],:,:) = rsq{3};
% rsq{3} = temp;
% temp = NaN(8,14,10);
% temp([1:3,5,7,8],:,:) = rsqout{3};
% rsqout{3} = temp;
% 
% temp = NaN(8,14,10);
% temp([1:3,5:8],:,:) = rsq{4};
% rsq{4} = temp;
% temp = NaN(8,14,10);
% temp([1:3,5:8],:,:) = rsqout{4};
% rsqout{4} = temp;
% save('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\SubspaceComparison\self_compare.mat','rsq','rsqout');

load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\SubspaceComparison\self_compare.mat');
%% example communication subspace
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ExampleCommunicationSubspaces';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotExampleSubspace(data,'VIS',5,1)
PlotExampleSubspace(data,'VIS',6,1)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleSubspaces_IN',savedir,0); close all
PlotExampleSubspace(dataout,'VIS',5,1)
PlotExampleSubspace(dataout,'VIS',6,1)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleSubspaces_OUT',savedir,0); close all
PlotExamplePredictors(data)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExamplePredictors',savedir,0); close all


%plot for all motifs
PlotExampleSubspace(data,'VIS',1:14,1)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleSubspaces_Combined',savedir,0); close all
%% Across all 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspaceDimensionality';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotSubspaceDimensionality(data,dataout,'VIS',1)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SubspaceDimensionality',savedir,0); close all

PlotSubspaceDimensionality(data,data,[],1)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SubspaceDimensionality_all',savedir,0); close all

%% See differences in beta weights across dimensions
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\BetaWeightsAcrossDimensions';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotDifferencesInBetaAcrossDimensions(data,'VIS', [5,14])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'BetaWeightsByDim_IN',savedir,0); close all

%% multiplexing
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Mulitplexing';
if ~exist(savedir,'dir'); mkdir(savedir); end
Plot_SubspaceMultiplexing(data)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'BetaWeightsByDim_IN',savedir,0); close all

%% Loading the sharing data
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\SubspaceComparison';
% folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\SubspaceComparison_meansubtract';
beta = cell(1,6);
betaout = cell(1,6);
for i = 1:6
    [fn,~] = GrabFiles(['type2_shuf1_rec',num2str(i),'.mat'],0,{folder}); 
    beta{i} = cellfun(@(x) load(x,'rsq'),fn);
    [fn,~] = GrabFiles(['type3_shuf1_rec',num2str(i),'.mat'],0,{folder}); 
    betaout{i} = cellfun(@(x) load(x,'rsq'),fn);
end
% adjust for the fact no RSP and bfd in rec 3 and 4
for i = 1:size(beta{3},2)
   temp = NaN(8,14,14,10);
   temp([1:3,5,7,8],:,:,:) = beta{3}(i).rsq;
   beta{3}(i).rsq = temp;
   temp = NaN(8,14,14,10);
   temp([1:3,5,7,8],:,:,:) = betaout{3}(i).rsq;
   betaout{3}(i).rsq = temp;
end
for i = 1:size(beta{4},2)
   temp = NaN(8,14,14,10);
   temp([1:3,5:8],:,:,:) = beta{4}(i).rsq;
   beta{4}(i).rsq = temp;
   temp = NaN(8,14,14,10);
   temp([1:3,5:8],:,:,:) = betaout{4}(i).rsq;
   betaout{4}(i).rsq = temp;
end
%%
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspaceSharing';
if ~exist(savedir,'dir'); mkdir(savedir); end
multicompresults = PlotSubspaceSharing(data,beta,'VIS',[5,14],1,data,rsq);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SubspaceSharingIN',savedir,0); close all


% multicompresults = PlotSubspaceSharing(dataout,betaout,'VIS',[5,14],1,data,rsqout);
% saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SubspaceSharingOut',savedir,0); close all

%% What do these dimensions mean? 

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspaceOrganization';
if ~exist(savedir,'dir'); mkdir(savedir); end
Plot_SubspaceOrganization(data); 
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'INorg',savedir,0); close all


%% Done 

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ModellingSubstractVersusDivisiveNormalization';
if ~exist(savedir,'dir'); mkdir(savedir); end
%how might different neurons be influenced by the normalization?
nval = 1000;
fr = linspace(0,5,nval);
a = NaN(nval,nval);
b = NaN(nval,nval);
for i = 1:nval
    for j = 1:nval
       a(i,j) = fr(i)/(fr(j)+1); 
       b(i,j) = fr(i)-fr(j); 
    end
end
figure; hold on; 
imagesc(a,[0 3]); colorbar; colormap magma
set(gca,'XScale','log','YScale','log')
tickval = get(gca,'xtick');
set(gca,'xtick',tickval,'xticklabel',round(fr(tickval),2));
set(gca,'ytick',tickval,'yticklabel',round(fr(tickval),2));
ylabel('Trial FR (spikes/sec)'); 
xlabel('Baseline FR (spikes/sec)'); 
title({'Divisive normalization increases','impact of sparse an selective neurons'},'fontweight','normal');
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'NormImpact',savedir,0); close all
%plot FR mean baseline vs FR mean post baseline versus beta weights 

%% Do they follow a power law? 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\PowerLawDimensionality';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotPowerLawDimensionality(data,1) 
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'INorg',savedir,0); close all

%% Compare the In and Out cortical layers

cur_a = 6;
yy = dataout;
cur_rec = 1;
[~,~,~,EphysPath] = LoadDataDirectories(cur_rec);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,yy{cur_rec}(1).st_depth);
neu_area = cat(2,neu_area{:});

temp = cat(1,neu_area(:).detailed_label);
[area_name, area_label] = ParseByArea(temp,neu_area,'general');
[area_name, area_label] = CleanUpAreas(area_name, area_label, 10); 

%load beta weights for each dimensions
b=[];
for cur_d = 1:10
%     [beta,area_label] = LoadVariable(yy,'rrr_V',area_label{cur_a},cur_d);
    [beta,area_label] = LoadVariable(yy,'rrr_beta',area_label{cur_a},cur_d);
    b(:,:,cur_d) = squeeze(beta(cur_rec,:,:));
end

%for each general area, plot the beta weights within that area (across dimensions)
peraarea = cellfun(@(x) sum(squeeze(nanmean(b(:,strcmp(area_name{cur_a},x),1:10),2)),2),unique(area_name{cur_a}),'UniformOutput',0);

x = cat(2,peraarea{:});
% y  = x;
figure; hold on
boxplot(x,'notch','on');
set(gca,'XTickLabel',unique(area_name{cur_a}),'XTickLabelRotation',45)
title('In subspace');
[p,~,stats] = kruskalwallis(x',[],'off');
title(sprintf('p=%0.3f',p),'fontweight','normal')




























