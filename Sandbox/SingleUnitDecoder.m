function [auc_vals,area_label,auc_vals_null] = SingleUnitDecoder(EphysPath,motif_fits,motifs,win, type)
%Camden - timeless
%get the accuracy of individual units at decoding between two motifs
%uses normalized fr during 1 second window after motif onset. 

if nargin <4; win = [1 15]; end %post onset during with which to average over
if nargin <5; type = 'peak'; end %avg or peak for the type of signal to use to decode the two stimuli


% EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
% [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
% motifs = [9,13];


%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%parse onsets
[~,trig_st_a] = ParseByOnset([],st_norm,motif_onset,win,motifs(1));
[~,trig_st_b] = ParseByOnset([],st_norm,motif_onset,win,motifs(2));

%merge across probes 
trig_st_a = cat(1,trig_st_a{:});
trig_st_b = cat(1,trig_st_b{:});

%get the summative activity
switch type
    case 'avg'
        trig_st_a = squeeze(nanmean(trig_st_a,2));
        trig_st_b = squeeze(nanmean(trig_st_b,2));
    case 'peak'
        trig_st_a = squeeze(nanmax(trig_st_a,[],2));
        trig_st_b = squeeze(nanmax(trig_st_b,[],2));                
    otherwise
        error('unknown activity type');
end

%subsample to balance
rng('default');
n = min(size(trig_st_a,2),size(trig_st_b,2));
trig_st_a = trig_st_a(:,randperm(size(trig_st_a,2),n));
trig_st_b = trig_st_b(:,randperm(size(trig_st_b,2),n));

predictors = cat(2,trig_st_a,trig_st_b)';
response = cat(1,ones(size(trig_st_a,2),1),2*ones(size(trig_st_b,2),1));
cvp = cvpartition(response, 'Holdout', 0.30);

%get decoding accuracy and auc for each neuron
auc_full = NaN(size(predictors,2),1);
for i = 1:size(predictors,2)
    [~, Observed, ~, ~] = SVMClassifier_Binary([predictors(:,i),response],cvp,...
        'nshuf',0,'pca',0,'solver',1,'kernel','linear','numkfold',10,'featureselect','none',...
        'optimize',0,'optimize_maxiter',100,'holdout',0.3);    
    auc_full(i) = Observed.AUC;
end

%parse the auc
[auc_vals, area_label] = ParseByArea(auc_full,neu_area,'parent');

%% get the null, permuted distribution of betas per neuron/area
rng('default'); 
auc_vals_null = cell(numel(auc_vals),1000);
for cur_perm = 1:1000
    cur_perm
    response = cat(1,ones(size(trig_st_a,2),1),2*ones(size(trig_st_b,2),1));
    response = response(randperm(numel(response),numel(response)));
    cvp = cvpartition(response, 'Holdout', 0.30);
    %get decoding accuracy and auc for each neuron
    auc_full = NaN(size(predictors,2),1);
    for i = 1:size(predictors,2)
        [~, Observed, ~, ~] = SVMClassifier_Binary([predictors(:,i),response],cvp,...
            'nshuf',0,'pca',0,'solver',1,'kernel','linear','numkfold',10,'featureselect','none',...
            'optimize',0,'optimize_maxiter',100,'holdout',0.3);    
        auc_full(i) = Observed.AUC;
    end
    %parse the betas
    [auc_vals_null(:,cur_perm), ~]  = ParseByArea(auc_full,neu_area,'parent');
end



end %function end


% 
% auc_full = NaN(size(predictors,2),1);
% for i = 1:size(predictors,2)
%     [~, Observed, ~, ~] = SVMClassifier_Binary([predictors(:,i),response],[],...
%         'nshuf',0,'pca',0,'solver',1,'kernel','linear','numkfold',10,'featureselect','none',...
%         'optimize',0,'optimize_maxiter',100,'holdout',0.3);    
%     auc_full(i) = Observed.AUC;
% end
% %%
% figure; hold on;
% c = getColorPalet(50);
% imagesc(auc_full',[0.5 1]); colormap(gca,flipud(gray))
% set(gca,'xlim',[0 numel(auc_full)],'XTick','','YTick','')
% PlotProbeAnatomy(gca, neu_area, 0.5, 'parent',0,0,c(randperm(size(c,1),size(c,1)),:)); 
















