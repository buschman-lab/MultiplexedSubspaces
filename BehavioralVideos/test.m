opts.detrend = 0 ; 
opts.fps = 60; 
opts.method = 'movingavg';

vid_name = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\RestingStateHemo\Mouse432_10_17_2019\Cam_0_20191017-171729.avi';
wf_name = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Long_Recordings\432-10-17-2019_1Fitted_block_hemoflag0_1.mat';

%load motif weightings
H = load(wf_name,'H');
H = H.H; 
num_frames = size(H,2);

%load and parse behavioral video
[avg_trace, data] = ParseVideos(vid_name);
data(1) = []; %remove timing signal since not used any more

%manual theshold for timing trace
[onset, offset] = ManualThreshold(avg_trace{1});

%Sliding Window DFF on the behavioral videos
w = opts.fps*30; 
dff = cellfun(@(x) makeDFF(x,opts,w), data, 'UniformOutput',0);
dff = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), dff, 'UniformOutput',0);

%PCA to get independent spatial modes of videos. pixels are variables
[spatial_modes, temporal_weights] = cellfun(@(x) pca(x','NumComponents',10),dff,'UniformOutput',0);

%get the average dff to compute just gross motor activity
avg_dff = cellfun(@(x) squeeze(nanmean(x,1:2)), dff,'UniformOutput',0);

%Trim 
avg_dff = cellfun(@(x) x(onset:offset),avg_dff,'UniformOutput',0);
temporal_weights = cellfun(@(x) x(onset:offset,:),temporal_weights,'UniformOutput',0);

%Downsample with linear interpolation to match length of imaging data
avg_dff_ds = cellfun(@(x) interp1(1:numel(x), x, linspace(1,numel(x),num_frames),'linear'),avg_dff,'UniformOutput',0);
temporal_weights_ds = cellfun(@(x) interp1(1:numel(x), x, linspace(1,numel(x),num_frames),'linear'),num2cell(temp,1),'UniformOutput',0);
temporal_weights_ds = cat(1,temporal_weights_ds{:});

%motif-triggered behavioral measures
[h_thresh, h_bin, h_onset] = ThresholdMatrix(H,3);


%% X Correlation between motifs and events
H_cent = H-nanmean(H,2);
%%
for j = 1:10
temp = [];

for i = 1:14
temp(i,:) = xcorr(H_cent(i,:)',temporal_weights_ds(j,:)',130);
end
figure; 
plot(temp')
end

%% Things to look at 
%1 gross motor activity
%2 specific PCA traces
%3 behavioral motifs

for i = 1:14
    idx1 = h_onset{i}<=1000;
    idx2 = h_onset{i}>=num_frames-1000;
    h_onset{i}(idx2==1)=[];
    h_onset{i}(idx1==1)=[];
    
end


%%
array_mean = [];
array_std = [];
temp = temporal_weights_ds(1,:);
val = 130;
for i = 1:numel(h_onset)
    figure();
    test = arrayfun(@(x) temp(x-val:x+val), h_onset{i},'UniformOutput',0);
    test = cat(1,test{:});
    array_mean(i,:) = nanmean(test,1);
    array_std(i,:) = nanstd(test,[],1)/sqrt(size(test,2));
    shadedErrorBar(1:val*2+1,array_mean(i,:),array_std(i,:))
   
   
end




%%YOU MAY WANT TO DIVIDE BY MOTOR ACTIVITY? BECAUSE THAT WILL IMPACT
%%THRESHODING

%correlation

%H spike triggered plot



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
COUNT = 1;
for a = [0,1]
    for b = [0,1]
        for c = [0,1]
            for d = [0.1, 0.01]
                for e = [50, 200]
                    params{COUNT} = [a,b,c,d,e];
                    COUNT = COUNT+1;
                end
            end
        end
    end
end

for cur_run = 1:numel(params)
cd('C:\Users\Camden\Documents\Work\OneDrive\Buschman Lab\Scratch Data\New Folder')
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

umap_py.min_dist = 0.9;
umap_py.n_neighbors = 50;
manifold_vals = umap_py.fit(projections(14400:2:floor(size(projections,1)/2),:));
fprintf('\n dimensionality reduction of %d timepoints and %d dimesions took %d minutes\n',...
    size(projections,1),size(projections,2),round(toc/60,0));
figure; 
plot(manifold_vals(:,1),manifold_vals(:,2),'linestyle','none','marker','.');
figure;
plot3(1:size(manifold_vals,1),manifold_vals(:,1),manifold_vals(:,2),'linestyle','none','marker','.','markersize',2,'color','k');

%% cluster the values with phenograph
[ovr_comm, ovr_Q] = PhenographSimple(manifold_vals, 'knn', bp.pheno_knn);
figure; plot(ovr_comm(:,1))

selftransform = (sum(diff(ovr_comm(:,1))==0))/numel(ovr_comm(:,1));
 
save('rundata.mat','ovr_comm','manifold_vals','params','selftransform','-v7.3')
close all;

end


%%

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

%motif-triggered behavioral measures
[h_thresh, h_bin, h_onset] = ThresholdMatrix(H_weight,3);

%% Resample the low dimensional space to match the frequency of the motifs
manifold_vals_resampled = NaN(num_frames,size(manifold_vals,2));
for i = 1:size(manifold_vals,2)
    temp = manifold_vals(:,i);
    manifold_vals_resampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');
end

%% for each motif, get the time bins around it

temp = h_onset{4};
val = 2*13;
temp(temp<=val)=[];
temp(temp>=(num_frames-val))=[];
dur = 2*13;
traj = arrayfun(@(x) manifold_vals_resampled(x-dur:x+dur,:), temp, 'UniformOutput',0);
traj_x = cellfun(@(x) x(:,1)', traj, 'UniformOutput',0);
traj_x = cat(1,traj_x{:});
traj_y = cellfun(@(x) x(:,2)', traj, 'UniformOutput',0);
traj_y = cat(1,traj_y{:});

figure; hold on; 
plot(manifold_vals(:,1),manifold_vals(:,2),'linestyle','none','marker','.');

for i = 1:numel(temp)
plot(traj_x(i,:),traj_y(i,:),'r','linewidth',2);
pause();
end





%phenograph to cluster these trajectories per motif


%get the number of behavioral states for each motif 


%%%%%%% SCRATCH CODE

% % % % % % %optional zscore each trace
% % % % % % if bp.zscore 
% % % % % %    projections = zscore(projections,1);
% % % % % % end
% % % % % % 
% % % % % % %wavelet transform
% % % % % % fprintf('\nCalculating Wavelet Transform\n');
% % % % % % [amps,f] = findWavelets(projections,bp); 
% % % % % % 
% % % % % % %30 second Sliding Window DFF on the behavioral videos
% % % % % % w = bp.fps*bp.window; 
% % % % % % %reshape to 3d to match dff function
% % % % % % facecam_mean_trim = cellfun(@(x) reshape(x,1,1,size(x,1)), facecam_mean_trim,'UniformOutput',0);
% % % % % % facecam_dff = cellfun(@(x) squeeze(makeDFF(x,bp,w)), facecam_mean_trim, 'UniformOutput',0);

%%
%gscatter(manifold_vals(:,1),manifold_vals(:,2),ovr_comm(:,1));
%make movie of the trajectory over time
temp = 1;
col = getColorPalet(numel(unique(ovr_comm(:,temp))));
figure; hold on; 
plot(manifold_vals(:,1),manifold_vals(:,2),'linestyle','none','marker','.','color',[0.8 0.8 0.8]);
for i = 1:size(manifold_vals,1)
    scatter(manifold_vals(i,1),manifold_vals(i,2),20,'r','.')
    drawnow();
end

















