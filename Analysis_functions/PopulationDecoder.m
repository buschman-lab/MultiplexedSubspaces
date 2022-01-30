function [beta,area_label,auc_val,acc_val,beta_null,auc_val_null,acc_val_null] = PopulationDecoder(EphysPath,motif_fits,motifs,win,type,n_perm)
%Camden - timeless
%get the accuracy across all recorded neurons and then determin each
%neurons beta weight%uses normalized fr during 1 second window after motif onset. 
%Input 'motifs' is a # pairings x 2 for motifs to compare

if nargin <4; win = [1 15]; end %post onset during with which to average over
if nargin <5; type = 'peak'; end %avg or peak for the type of signal to use to decode the two stimuli

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct

%loop through each pair of motifs
beta = cell(size(motifs,1),1);
beta_null = cell(size(motifs,1),1);
area_label = cell(size(motifs,1),1);
auc_val = NaN(size(motifs,1),1);
acc_val = NaN(size(motifs,1),1);
auc_val_null = NaN(size(motifs,1),n_perm);
acc_val_null = NaN(size(motifs,1),n_perm);
for cur_pair = 1:size(motifs,1)
    %parse onsets
    [~,trig_st_a] = ParseByOnset([],st_norm,motif_onset,win,motifs(cur_pair,1));
    [~,trig_st_b] = ParseByOnset([],st_norm,motif_onset,win,motifs(cur_pair,2));

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

    %% decode across neurons
    predictors = cat(2,trig_st_a,trig_st_b)';
    response = cat(1,ones(size(trig_st_a,2),1),2*ones(size(trig_st_b,2),1));
    cvp = cvpartition(response, 'Holdout', 0.30);
    [~, Observed, ~, ~] = SVMClassifier_Binary([predictors,response],cvp,...
        'nshuf',0,'pca',0,'solver',1,'kernel','linear','numkfold',4,'featureselect','none',...
        'optimize',0,'optimize_maxiter',100);

    auc_val(cur_pair) = Observed.AUC;
    acc_val(cur_pair) = Observed.Accuracy;

    %parse the betas
    [beta{cur_pair}, area_label{cur_pair}] = ParseByArea(Observed.Classifier.Beta,neu_area,'parent');

    %% get the null, permuted distribution of betas per neuron/area
    if n_perm >0
        rng('default'); 
        beta_null = cell(numel(beta{cur_pair}),n_perm);
        stats_null = [];
        for i = 1:n_perm
            response = cat(1,ones(size(trig_st_a,2),1),2*ones(size(trig_st_b,2),1));
            response = response(randperm(numel(response),numel(response)));
            cvp = cvpartition(response, 'Holdout', 0.30);
            [~, Observed, ~, ~] = SVMClassifier_Binary([predictors,response],cvp,...
                'nshuf',0,'pca',0,'solver',1,'kernel','linear','numkfold',4,'featureselect','none',...
                'optimize',0,'optimize_maxiter',100);
            %parse the betas
            [beta_null{cur_pair}(:,i), ~] = ParseByArea(Observed.Classifier.Beta,neu_area,'parent');    
            auc_val_null(cur_pair,i) = Observed.AUC;
            acc_val_null(cur_pair,i) = Observed.Accuracy;
        end
    else
        stats_null = [];
        beta_null = [];
    end
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
















