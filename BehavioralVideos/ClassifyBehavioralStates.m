function ClassifyBehavioralStates()
%Camden MacDowell - timeless


%Load the processed weights from each animal (see BehavioralState_Analysis)
weight_431 = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse431_10_17_2019\processed_weights.mat');
weight_431 = weight_431.weight;
weight_432 = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse432_10_17_2019\processed_weights.mat');
weight_432 = weight_432.weight;

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

%now classify behavioral state by the motif values
val = ([5, 10, 15, 20, 25, 30])/100;
auc = [];
for i = 1:numel(val)
rng('default');
inactive =  cat(1,weight_norm_state{:,1})'; %state 1
inactive(any(isnan(inactive),2),:)=[]; %remove NaNs
active =  cat(1,weight_norm_state{:,2})'; %state 2
active(any(isnan(active),2),:)=[]; %remove NaNs
num_sample = min(size(active,1),size(inactive,1));
%get equal samples from both classes
active = active(randperm(size(active,1),num_sample),:);
inactive = inactive(randperm(size(inactive,1),num_sample),:);
data = [cat(1,inactive,active), cat(1,ones(size(inactive,1),1),2*ones(size(active,1),1))];
[~, Observed, ~, TrainAUC] = SVMClassifier_Binary(data,[],'holdout',0.1,...
    'nshuf',0,'pca',0,'solver',1,'kernel','rbf','numkfold',15,'featureselect','none',...
    'optimize',1,'optimize_maxiter',75);
auc(i) = Observed.AUC;
end


  
















