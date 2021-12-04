function CompareMotifDiscoveryParameters_Figures()
%Camden MacDowell - timeless
%Goes with 'CompareMotifDiscoveryParameters.m'
%right now, only designed to compare length

%loop through the folders in this directory
fold_list = dir('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\MotifDiscoveryParams');
fold_list = fold_list(3:end); 

file_list = cell(1,numel(fold_list));
for i = 1:numel(fold_list)
   file_list{i} = GrabFiles('\w*chunk\w*.mat',0,{[fold_list(i).folder,filesep, fold_list(i).name]});
end

%% compare statistics 
pev = cell(1,numel(file_list));
n_motifs = cell(1,numel(file_list));
cost = cell(1,numel(file_list));
pev_test = cell(1,numel(file_list));
n_motifs_test = cell(1,numel(file_list));
cost_test = cell(1,numel(file_list));
basis_n = NaN(1,numel(file_list));
basis_q = NaN(1,numel(file_list));
for i = 1:numel(file_list)
   temp = cellfun(@(x) load(x,'stats_train'),file_list{i},'UniformOutput',0); 
   train_stats = cellfun(@(x) x.stats_train, temp,'UniformOutput',0); 
   pev{i} = cellfun(@(x) x.pev, train_stats,'UniformOutput',1); 
   n_motifs{i} = cellfun(@(x) x.n_motifs, train_stats,'UniformOutput',1); 
   cost{i} = cellfun(@(x) x.cost, train_stats,'UniformOutput',1);    
   temp = cellfun(@(x) load(x,'stats_test'),file_list{i},'UniformOutput',0); 
   test_stats = cellfun(@(x) x.stats_test, temp,'UniformOutput',0); 
   pev_test{i} = cellfun(@(x) x.pev, test_stats,'UniformOutput',1); 
   n_motifs_test{i} = cellfun(@(x) x.n_motifs, test_stats,'UniformOutput',1); 
   cost_test{i} = cellfun(@(x) x.cost, test_stats,'UniformOutput',1);    
   %basis motif stats
   if exist([fold_list(i).folder,filesep, fold_list(i).name filesep,'_basis_motifs.mat'])
       temp = load([fold_list(i).folder,filesep, fold_list(i).name filesep,'_basis_motifs.mat'],'ovr_q','cluster_idx');
       basis_n(i) = numel(unique(temp.cluster_idx));
       basis_q(i) = temp.ovr_q;
   end
end

%get the order
motif_L = NaN(1,numel(fold_list));
for i = 1:numel(fold_list)
   motif_L(i) = str2double(erase(fold_list(i).name,'L'));
end

% recorder
[~, idx] = sort(motif_L,'ascend');

%% figures
figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, n_motifs(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), n_motifs(idx),'UniformOutput',1),'lineprops',{'color','r'})

figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, cost(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), cost(idx),'UniformOutput',1),'lineprops',{'color','r','marker','o'})

figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, pev(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), pev(idx),'UniformOutput',1),'lineprops',{'color','r','marker','o'})

figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, pev_test(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), pev_test(idx),'UniformOutput',1),'lineprops',{'color','r','marker','o'})

figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, n_motifs_test(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), n_motifs_test(idx),'UniformOutput',1),'lineprops',{'color','r'})

figure; hold on; 
shadedErrorBar(motif_L(idx),cellfun(@nanmean, cost_test(idx),'UniformOutput',1),cellfun(@(x) sem(x,2), cost_test(idx),'UniformOutput',1),'lineprops',{'color','r','marker','o'})

figure; hold on; 
plot(motif_L(idx),basis_n(idx),'color','r','marker','o')

figure; hold on; 
plot(motif_L(idx),basis_q(idx),'color','r','marker','o')

%% motif duration
for i = 1:numel(file_list)
   PlotMotifDuration(file_list{i})
end


%maybe plot the cost as well? and PEV and generalized PEV. And equity
%between the the number of motifs in the testing

