function ClassifyBehavioralStates()
%Camden MacDowell - timeless
%Load the processed weights from each animal (see BehavioralState_Analysis)
if ispc
    weight_431 = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse431_10_17_2019\processed_weights.mat');    
    weight_432 = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse432_10_17_2019\processed_weights.mat');
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    weight_431 = load('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/Mouse431_10_17_2019/processed_weights.mat');    
    weight_432 = load('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/Mouse432_10_17_2019/processed_weights.mat');    
end
weight_432 = weight_432.weight;
weight_431 = weight_431.weight;

weight =  weight_431;
for i = 1:size(weight,1)
   for j = 1:size(weight,2)
       weight{i,j} = [weight_431{i,j},weight_432{i,j}];
   end       
end
%first normalize by the column to get relative motif weighting vs baseline
avg_weight = [nanmean(([weight{:,1}])),nanmean(([weight{:,2}]))];
weight_norm = weight;
for i = 1:size(weight,1)
   for j = 1:size(weight,2)
       weight_norm{i,j} = (weight{i,j})/avg_weight(j);
   end       
end
%first normalize by the column to get relative motif weighting across the two states
weight_norm_state = weight_norm;
for i = 1:size(weight_norm,1)
    temp = cat(1,weight_norm_state{i,:})'; 
    for j = 1:size(weight_norm,2)
        weight_norm_state{i,j} = weight_norm_state{i,j}/nanmean(temp(:));
    end         
end

% organize data
rng('default'); %for reproducibility
inactive =  cat(1,weight_norm_state{:,1})'; %state 1
inactive(any(isnan(inactive),2),:)=[]; %remove NaNs
active =  cat(1,weight_norm_state{:,2})'; %state 2
active(any(isnan(active),2),:)=[]; %remove NaNs
num_sample = min(size(active,1),size(inactive,1));
%get equal samples from both classes
active = active(randperm(size(active,1),num_sample),:);
inactive = inactive(randperm(size(inactive,1),num_sample),:);

%classify by the full space
fprintf('\n\tClassifying by the full space\n');
rng('default'); %for reproducibility
data = [cat(1,inactive,active), cat(1,ones(size(inactive,1),1),2*ones(size(active,1),1))];
[~, Observed, ~, ~] = SVMClassifier_Binary(data,[],'holdout',0.2,...
    'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',10,'featureselect','none',...
    'optimize',1,'optimize_maxiter',100);
auc_full = Observed.AUC;

%leave-one-out
auc_leave = NaN(1,size(inactive,2));
for cur_motif = 1:size(inactive,2)
    rng('default'); %for reproducibility
    idx = ones(1,size(inactive,2));
    idx(cur_motif)=0;
    fprintf('\n\tClassifying by the full space\n');
    data = [cat(1,inactive(:,idx==1),active(:,idx==1)), cat(1,ones(size(inactive,1),1),2*ones(size(active,1),1))];
    [~, Observed, ~, ~] = SVMClassifier_Binary(data,[],'holdout',0.2,...
        'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',10,'featureselect','none',...
        'optimize',1,'optimize_maxiter',100);
    auc_leave(cur_motif) = Observed.AUC;
end

%classify by the full space again with multiple randomizations to get a distirbution of success rates
fprintf('\n\t iterating classification\n');
rng('default'); %for reproducibility
auc_full_iter = NaN(1,100);
for i = 1:100
    data = [cat(1,inactive,active), cat(1,ones(size(inactive,1),1),2*ones(size(active,1),1))];
    [~, Observed, ~, ~] = SVMClassifier_Binary(data,[],'holdout',0.2,...
        'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',10,'featureselect','none',...
        'optimize',1,'optimize_maxiter',100);
    auc_full_iter(i) = Observed.AUC;
end

%get the contibutions of each motif to classification accuracy
if ~ispc
    if ~exist('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/LeaveOneOutClassification/')
        mkdir('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/LeaveOneOutClassification/');
    end
    save('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/LeaveOneOutClassification/LeaveOneOut.mat','auc_leave','auc_full','auc_full_iter')
end

end %function

% %now classify behavioral state by the motif values
% param_mat = combvec((0.1:0.05:0.4),(3:2:20));
% 
% rng('default');
% inactive =  cat(1,weight_norm_state{:,1})'; %state 1
% inactive(any(isnan(inactive),2),:)=[]; %remove NaNs
% active =  cat(1,weight_norm_state{:,2})'; %state 2
% active(any(isnan(active),2),:)=[]; %remove NaNs
% num_sample = min(size(active,1),size(inactive,1));
% %get equal samples from both classes
% active = active(randperm(size(active,1),num_sample),:);
% inactive = inactive(randperm(size(inactive,1),num_sample),:);data = [cat(1,inactive,active), cat(1,ones(size(inactive,1),1),2*ones(size(active,1),1))];
% [~, Observed, ~, trained] = SVMClassifier_Binary(data,[],'holdout',param_mat(1,block),...
%     'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',param_mat(2,block),'featureselect','none',...
%     'optimize',0,'optimize_maxiter',100);
% 
% if ~ispc
%     if ~exists('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/linearclassification/')
%         mkdir('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/linearclassification/');
%     end
%     save(['/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData_MesomappingManuscript_5_2019/DeepLabCut_BehavioralState_Analysis/linearclassification/',sprintf('Block%d.mat',block)],'Observed','trained')
% end
%     
%
% %%parameter sweep plotting function
% files = GrabFiles('Block',0,{pwd});
% auc = NaN(1,numel(files));
% for i = 1:numel(files)
%     temp = load(files{i});
%     auc(i) = temp.Observed.AUC;
% end


  
















