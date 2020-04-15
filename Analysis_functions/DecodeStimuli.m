function DecodeStimuli(rep,varargin)
%Camden MacDowell - timeless
%load the maximum sensory response per trial split by motif and stimulus 

if ispc
    data = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\SensoryRepitoires\PeakWeightMotifPerStim2Seconds.mat','PeakWeight');
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    data = load('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/PeakWeightMotifPerStim2Seconds.mat','PeakWeight');
%     data = load('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/PeakWeightMotifPerStim1Seconds.mat','PeakWeight');
end
data = data.PeakWeight;

opts.optimize = 0; 
opts.holdout = 0.2;
opts.removenan = 0; %remove trials that failed to evoked all motifs. Def=0
opts = ParseOptionalInputs(opts,varargin);

%randomize rng. 
rng('shuffle'); 

nM = numel(data);

%% Decoding accuracy of full model and LOO
% Balanced by motif AND stimuli

if opts.removenan
    %get the minimum samples across all motifs
    num_sample = NaN(nM,1);
    for cur_motif = 1:nM 
        x = data(cur_motif).motif{1};    
        y = data(cur_motif).motif{2};
        x(isnan(x))=[];
        y(isnan(y))=[];
        num_sample(cur_motif) = min(numel(y),numel(x));    
    end %data compiling
    num_sample = min(num_sample);
end

%Compile the data
full_model = [];
for cur_motif = 1:nM 
    x = data(cur_motif).motif{1};
    y = data(cur_motif).motif{2};
    
    if opts.removenan
        x(isnan(x))=[];
        y(isnan(y))=[];
        x = x(randperm(numel(x),num_sample));
        y = y(randperm(numel(y),num_sample)); 
    else
        x(isnan(x))=0;
        y(isnan(y))=0;
    end    
    full_model(:,cur_motif) = cat(1,x,y);
end %data compiling
response = cat(1,ones(size(x)),2*ones(size(y))); 

%correlation the trial-by-trial responses to motif 1 and 4; 
cvp = cvpartition(response, 'Holdout', opts.holdout);

%Classify full model
[~, Observed, ~, ~] = SVMClassifier_Binary([full_model,response],cvp,...
    'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',5,'featureselect','none',...
    'optimize',opts.optimize,'optimize_maxiter',50);
auc_full = Observed.AUC;

%Classify LOO
auc_loo = NaN(1,nM);
for cur_motif = 1:nM
    loo_model = full_model;
    loo_model(:,cur_motif) = [];
    [~, Observed, ~, ~] = SVMClassifier_Binary([loo_model,response],cvp,...
        'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',5,'featureselect','none',...
        'optimize',opts.optimize,'optimize_maxiter',50);
    auc_loo(cur_motif) = Observed.AUC;
end

% Decoding accuracy of each motif.
auc_per_motif = NaN(1,nM);
for cur_motif = 1:nM        
    x = data(cur_motif).motif{1};    
    y = data(cur_motif).motif{2};        
    if opts.removenan
        x(isnan(x))=[];
        y(isnan(y))=[];
        num_sample = min(numel(y),numel(x));
        x = x(randperm(numel(x),num_sample));
        y = y(randperm(numel(y),num_sample)); 
    else
        x(isnan(x))=0;
        y(isnan(y))=0;
    end     
    response = cat(1,ones(size(x)),2*ones(size(y)));    
    
    %classify    
    cvp = cvpartition(response, 'Holdout', opts.holdout);
    [~, Observed, ~, ~] = SVMClassifier_Binary([cat(1,x,y),response],cvp,...
        'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',5,'featureselect','none',...
        'optimize',opts.optimize,'optimize_maxiter',50);
    
    auc_per_motif(cur_motif) = Observed.AUC;
end %motif decoding accuracy


%get the contibutions of each motif to classification accuracy
if ~ispc
    if ~exist('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/SensoryDecoderByMotif/','dir')
        mkdir('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/SensoryDecoderByMotif/');
    end
%     save(['/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/SensoryDecoderByMotif/',...
%         sprintf('1sec_removenan_%d_optimized_%d_holdout_%0.2g_replication_%g.mat',opts.removenan,opts.optimize,opts.holdout,rep)],'auc_per_motif','auc_full','auc_loo')
    save(['/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/SensoryRepitoires/SensoryDecoderByMotif/',...
        sprintf('removenan_%d_optimized_%d_holdout_%0.2g_replication_%g.mat',opts.removenan,opts.optimize,opts.holdout,rep)],'auc_per_motif','auc_full','auc_loo','opts')
end

end %function


% rng('default')
% tact = full_model(response==1,[1,4]);
% vis = full_model(response==2,[1,4]);
% for iter = 1:1000 %permuataion 
%     %get baseline    
%     base = randsample(size(vis,1),size(vis,1),'true') %visual is stim 2
%     rho_base = corr(vis(base,1),vis(base,2));
%     
%     
%     test = randsample( %tactile is stim 1
%     










