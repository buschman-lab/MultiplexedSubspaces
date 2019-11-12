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




















