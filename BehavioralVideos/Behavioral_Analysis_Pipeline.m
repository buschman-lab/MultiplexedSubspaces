%add paths
addpath(genpath('C:\Users\Camden\Documents\Github\MotionMapper'));
addpath(genpath('C:\Users\Camden\Documents\Github\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\Camden\Documents\Github\umap'));

%set filepaths
fn_path = 'C:\Users\Camden\Desktop\Videos\Mouse9031_10_17\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_dlc = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse9031_10_17';
fn_widefield = '431-10-17-2019_1Fitted_block_hemoflag0_1.mat';

condapath = 'C:\ProgramData\Anaconda3\Scripts\conda.exe';
condaenv = 'work_umap';

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

%% Precessing Behavioral Videos 
%parse the facecam to get timing signal and behavioral features
[facecam_mean, ~] = ParseVideos([fn_path, fn_facecam],bp);

%save the roi image
saveCurFigs(gcf,'-svg','Chosen ROIs',fn_path,1);
close; 

%manual theshold for timing trace
[onset, offset] = ManualThreshold(facecam_mean{1});

%Trim to match start and stop of imaging 
facecam_mean_trim = cellfun(@(x) x(onset:offset),facecam_mean(2:3),'UniformOutput',0);

%% Processing DLC Analyzed Data
[~,~,raw_data] = xlsread([fn_path, fn_dlc]);
%parse the dlc data by desired bodyparts, fitler, and reference to tail
[projections, Perc_Filtered] = parse_dlc(raw_data,bp);

%Trim to match start and stop of imaging 
projections = projections(onset:offset,:);

if bp.include_facecam %add the other behavioral traces to the dlc data   
   projections = [facecam_mean_trim{:} projections];
end

if bp.derivative %Get the absolute derivative (this is usually not good at all)
    projections = abs(diff(projections,1));
end
 
if bp.zscore %optional zscore
   projections = zscore(projections,1);
end

%umap dimensionality reduction in python
conda.init(condapath);
conda.setenv('work_umap')
umap_py = umap; 
fprintf('\n dimensionality reduction with umap')

umap_py.min_dist = bp.umap_mind_dist;
umap_py.n_neighbors = bp.umap_n_neighbors;
% manifold_vals = umap_py.fit(projections(14400:2:floor(size(projections,1)/5),:));
manifold_vals = umap_py.fit(projections);
fprintf('\n dimensionality reduction of %d timepoints and %d dimesions took %d minutes\n',...
    size(projections,1),size(projections,2),round(toc/60,0));

figure; 
plot(manifold_vals(:,1),manifold_vals(:,2),'linestyle','none','marker','.','markersize',1,'color',[0.5 0.5 0.5]);

%save off
save('umap_data.mat','manifold_vals', 'bp', '-v7.3')

%% load motif H weightings (optionally scale by weightings?)
%load motif weightings
H = load([fn_path, fn_widefield],'H');
H = H.H; 
num_frames = size(H,2);

%Get the average weight of each factor
load('C:\Users\Camden\Desktop\Videos\Mouse9031_10_17\AverageDPs_1.mat')
weights = (nanmean(W_clust_smooth,1));
H_weight = NaN(size(H));
for i = 1:size(H,1)
    H_weight(i,:) = helper.reconstruct(weights(:,i,:),H(i,:));
end

%% Resample the low dimensional space to match the frequency of the motifs
manifold_vals_resampled = NaN(num_frames,size(manifold_vals,2));
for i = 1:size(manifold_vals,2)
    temp = manifold_vals(:,i);
    manifold_vals_resampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');
end

movecam = facecam_mean_trim{2};
facecam_movemet = NaN(num_frames,size(movecam,2));
for i = 1:size(movecam,2)
    temp = movecam(:,i);
    movecam_resampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');
end

%% cluster the values with phenograph
[ovr_comm_resampled, ovr_Q] = PhenographSimple(manifold_vals_resampled, 'knn', bp.pheno_knn);

%save off
save('phenograph_data_downsampled.mat','ovr_comm_resampled', 'ovr_Q', '-v7.3')


%% for each motif, get the time bins around it
H_weight = H_weight(expvaridx,:);
[h_thresh, h_bin, h_onset] = ThresholdMatrix(H_weight,bp.min_std,bp.max_std);


motif_clust = {};
motif_dist_within={};
motif_dist_between={};
for cur_motif = 1:size(H_weight,1)
    temp = h_onset{cur_motif};
    temp(temp<=5)=[];
    temp(temp>=(num_frames-bp.trig_dur))=[];    
    
    %number of clusters moved in positional manifold 
    cluster_trajectory = arrayfun(@(x) ovr_comm_resampled(x-5:x+bp.trig_dur,bp.pheno_level), temp, 'UniformOutput',0);    
    unique_clusters_num = cellfun(@(x) numel(unique(x)), cluster_trajectory,'UniformOutput',0);
    
    %trajectory moved in positional manifold 
    manifold_trajectory = arrayfun(@(x) manifold_vals_resampled(x-5:x+bp.trig_dur,:), temp, 'UniformOutput',0); 
    manifold_distance = cellfun(@(x) sum(abs(diff(x,1)),'all'), manifold_trajectory, 'UniformOutput',0); 
    
    manifold_trajectory = cellfun(@(x) [[0;0] abs(diff(x,1))'], manifold_trajectory, 'UniformOutput',0); 
    manifold_distance_within = cellfun(@(x,y) sum(y(:,[0 diff(x)']==0),'all'), cluster_trajectory,manifold_trajectory,'UniformOutput',0);
    manifold_distance_between = cellfun(@(x,y) sum(y(:,[0 diff(x)']~=0),'all'), cluster_trajectory,manifold_trajectory,'UniformOutput',0);

    %parse data to plot
    motif_clust{cur_motif} = [unique_clusters_num{:}];   
    motif_dist_within{cur_motif} = [manifold_distance_within{:}];
    motif_dist_between{cur_motif} = [manifold_distance_between{:}];
end    

motif_clust = MakeCellsEqual(motif_clust,2,1);
motif_dist_within = MakeCellsEqual(motif_dist_within,2,1);
motif_dist_between = MakeCellsEqual(motif_dist_between,2,1);

motif_dist_anova_w = anova1(cat(1,motif_dist_within{:})',[],'off');
motif_dist_anova_b = anova1(cat(1,motif_dist_between{:})',[],'off');
motif_clust_anova = anova1(cat(1,motif_clust{:})',[],'off');

%% make a figure of motifs sorted by distance travelled within clusters
[~,idx] = sort(nanmean(cat(1,motif_dist_within{:}),2),'descend');
temp = cat(1,motif_dist_within{:});
temp = temp(idx,:);

col = getColorPalet(size(temp,1));
col = col(idx,:);
figure('position',[1092 252 500 586]); hold on
for i = 1:size(temp,1)
    y = temp(i,:);
    y(isnan(y))=[];
    y_sem = std(y)/sqrt(numel(y));
    errorbar(i,nanmean(y), y_sem,'color',col(i,:)','linewidth',2,'marker','.','markersize',20);    
end

%add pval
p = anova1(temp',[],'off');
AddSig(1,p,[1 size(temp,1) 42],1,0.15,1);
legend(num2str(idx),'Location','eastoutside')
set(gca,'XTick',[]);
ylabel({'Distance Traveled';'Within Clusters (AU)'})
xlabel({'Motif'});
ylim([38 42]);
xlim([0 15]);
setFigureDefaults();
set(gca,'position',[3 4 6 8.5]);


%% make a figure of motifs sorted by distance travelled between clusters
[~,idx] = sort(nanmean(cat(1,motif_dist_between{:}),2),'descend');
temp = cat(1,motif_dist_between{:});
temp = temp(idx,:);

col = getColorPalet(size(temp,1));
col = col(idx,:);
figure('position',[1092 252 500 586]); hold on
for i = 1:size(temp,1)
    y = temp(i,:);
    y(isnan(y))=[];
    y_sem = std(y)/sqrt(numel(y));
    errorbar(i,nanmean(y), y_sem,'color',col(i,:)','linewidth',2,'marker','.','markersize',20);    
end

%add pval
p = anova1(temp',[],'off');
AddSig(1,p,[1 size(temp,1) 150],1,4,1);
legend(num2str(idx),'Location','eastoutside')
set(gca,'XTick',[]);
ylabel({'Distance Traveled';'Between Clusters (AU)'})
xlabel({'Motif'});
ylim([80 150]);
xlim([0 15]);
setFigureDefaults();
set(gca,'position',[3 4 6 8.5]);


%% make a figure of motifs sorted by clusters travelled
[~,idx] = sort(nanmean(cat(1,motif_clust{:}),2),'descend');
temp = cat(1,motif_clust{:});
temp = temp(idx,:);

col = getColorPalet(size(temp,1));
col = col(idx,:);
figure('position',[1092 252 500 586]); hold on
for i = 1:size(temp,1)
    y = temp(i,:);
    y(isnan(y))=[];
    y_sem = std(y)/sqrt(numel(y));
    errorbar(i,nanmean(y), y_sem,'color',col(i,:)','linewidth',2,'marker','.','markersize',20);    
end

%add pval
p = anova1(temp',[],'off');
AddSig(1,p,[1 size(temp,1) 18],1,0.4,1);
legend(num2str(idx),'Location','eastoutside')
set(gca,'XTick',[]);
ylabel({'Number of Clusters'})
xlabel({'Motif'});
ylim([12 18]);
xlim([0 15]);
setFigureDefaults();
set(gca,'position',[3 4 6 8.5]);





%%
% gscatter(manifold_vals(:,1),manifold_vals(:,2),ovr_comm(:,1));
%make movie of the trajectory over time
temp = 1;
col = getColorPalet(numel(unique(ovr_comm_resampled(:,temp))));
figure; hold on; 
plot3((1:size(manifold_vals_resampled,1)),manifold_vals_resampled(:,1),manifold_vals_resampled(:,2),'linestyle','none','marker','.','color',[0.8 0.8 0.8]);
COUNT = 0;
% for i = 1:size(
%     COUNT = COUNT+1;
%     scatter3(x(i),manifold_vals_resampled(i,1),manifold_vals_resampled(i,2),20,col(ovr_comm_resampled(COUNT,temp),:),'.')
% end
% gscatter(manifold_vals_resampled(:,1),manifold_vals_resampled(:,2),ovr_comm_resampled(:,),[],[],0.01);
legend off