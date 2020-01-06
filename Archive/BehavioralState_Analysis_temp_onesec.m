%add paths
addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\'));
% % %set filepaths
fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse431_10_17_2019\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_bodycam = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
fn_dlc = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse431_10_17_2019';
fn_widefield = '431-10-17-2019_1Fitted_block_hemoflag0_1';
split_idx = {12,11.00,0.3}; %mouse 431

% % %set filepaths
% fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse432_10_17_2019\';
% fn_facecam = 'Cam_0_20191017-171729.avi';
% fn_bodycam = 'Cam_1_20191017-171729_Mouse432_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
% fn_dlc = 'Cam_1_20191017-171729_Mouse432_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
% fn_savebase = 'Mouse432_10_17_2019';
% fn_widefield = '432-10-17-2019_1Fitted_block_hemoflag0_1';
% split_idx = {11.58,6,0.3}; %mouse 432


%set filepaths
% fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse494_10_17_2019\';
% fn_facecam = 'Cam_0_20191017-184642.avi';
% fn_bodycam = 'Cam_1_20191017-184642_Mouse494_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
% fn_dlc = 'Cam_1_20191017-184642_Mouse494_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
% fn_savebase = 'Mouse494_10_17_2019';
% fn_widefield = '494-10-17-2019_1Fitted_block_hemoflag0_1';
% split_idx = {12.06,17.53,0.15};

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

H = load([fn_path, fn_widefield],'H');
num_frames = size(H.H,2);

savefigs = 0; 

%% Preprocessing Behavioral Videos 
%parse the facecam to get timing signal and behavioral features
if exist([fn_path fn_savebase 'roi_info.mat'],'file')
    roivals = load([fn_path fn_savebase 'roi_info.mat']);
    [facecam_mean, facecam_data, ~] = ParseVideos([fn_path, fn_facecam],bp,[],roivals.roi);
    onset = roivals.onset;
    offset = roivals.offset;
else
    [facecam_mean, facecam_data, roi] = ParseVideos([fn_path, fn_facecam],bp);
    %save the roi image
    axis off
    saveCurFigs(handles,'-svg','Chosen ROIs',fn_path,1);
    close all; 

    %manual theshold for timing trace
    [onset, offset] = ManualThreshold(facecam_mean{1});

    %save off the roi information 
    save([fn_path fn_savebase 'roi_info.mat'],'roi','onset','offset');
end

%get the average motion energy (absolute derivative) of the pixels
temp = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), facecam_data(2:end),'UniformOutput',0);
face_motion_energy = cellfun(@(x) mean(abs(diff(x,2)),1)', temp,'UniformOutput',0);

%% Processing DLC Analyzed Data
%load dlc traced data
[~,~,raw_data] = xlsread([fn_path, fn_dlc]);

%get the speed of the 4 limbs
[limbs, id, ~] = parse_dlc(raw_data,{'frontrightpawcenter','frontleftpawcenter','backrightpawcenter','backleftpawcenter'},[],bp.dlc_epsilon);
limb_speed = cellfun(@(x) [0; mean(abs(diff(limbs(:,strcmp(id,x)),1)),2)],{'frontrightpawcenter','frontleftpawcenter','backrightpawcenter','backleftpawcenter'},'UniformOutput',0); 
limb_speed = [limb_speed{:}];
limb_speed = mean(limb_speed,2);

%combined all features
% face_motion_energy{2} = -1 * face_motion_energy{2}+max(face_motion_energy{2}(:)); %may need to flip the whisker energy if high whisking actually blurs the camera and make low energy  
features = cat(2,face_motion_energy{:},limb_speed);
labels = {'nose motion energy','whisker motion energy','limb speed'};
labels_abbrev = {'NME','WME','LS'};

for i = 1:size(features,2)
    features(:,i) = convn(features(:,i),ones(130,1)/130,'same');
end

%Trim to match start and stop of imaging 
features = features(onset:offset,:);

% Downsample to match motif duration
features_downsampled = NaN(num_frames,size(features,2));
for i = 1:size(features_downsampled,2)
    temp = features(:,i);
    features_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear'); 
end

%store the mapping from original to downsampled 
x_query_ds = linspace(1,size(features,1),num_frames);

if bp.zscore %optional zscore
   features_downsampled = zscore(features_downsampled,1);
end

clear raw_data facecam_data W_clust_smooth
%% Parse states
num_features = size(features_downsampled,2); 
col = getColorPalet(num_features);

data = features_downsampled(42000:43300,:);
figure('units','centimeters','position',[1 1 15 25]); hold on;
ax =[];
x = (0:1300)/13;
for i = 1:size(data,2)
    ax(i) = subplot(size(data,2),1,i); 
    plot(x,data(:,i),'linewidth',2,'color',col(i,:))
    ylabel('Z-Score')
    if i == size(data,2)
        xlabel('Time (s)');
    else
        set(gca,'XColor','none') 
    end
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;       
end
for i = 1:size(data,2)
    set(ax(i),'position',[3 15-((i-1)*4) 10 3])
end

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Exampled_Traces',fn_path,1);
    close all;
end


figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
features_binned = NaN(size(features_downsampled));
x = linspace(0,60,size(features_downsampled,1));
for i = 1:num_features
    temp = NaN(size(features_downsampled(:,i)));
    temp(features_downsampled(:,i)<=split_idx{i})=0;
    temp(features_downsampled(:,i)>split_idx{i})=1;
    features_binned(:,i) = temp;    
    subplot(r,c,i); hold on;
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')  
    plot(x,temp,'color',[0.5 0.5 0.5],'linestyle','none','marker','.','markersize',5)
    set(gca,'ylim',[-0.25 1.25]);
    xlabel('Time (min)');
    ylabel('State');
    set(gca,'YTick',[0 1],'YTickLabels',{'low','high'})
    setFigureDefaults
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 6 6])
    
end

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','MetricsOverTime',fn_path,1);
    close all;
end


%distributions
figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
for i = 1:num_features
    subplot(r,c,i); hold on;
    %plot pdf 
    temp = features_downsampled(:,i);
    [f,xi] = ksdensity(temp); 
    plot(xi,f,'linewidth',2,'color',col(i,:)); 
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;
    ylabel('PDF')
    xlabel('Z-Score')
    yvals = get(gca,'ylim');
    line([split_idx{i},split_idx{i}],[0 max(get(gca,'ylim'))],'linestyle','--','linewidth',2,'color','k')
    set(gca,'ylim',yvals);
    if i ==4
        xlim([-1 4])
    end
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 3 3])
end

figure('position',[680   101   858   877]); hold on; 
[n,c] = numSubplot(size(features_downsampled,2),2);
x = linspace(0,60,size(features_downsampled,1));
for i = 1:size(features_downsampled,2)
    subplot(n,c,i); 
    plot(x,features_downsampled(:,i),'color',[0.5 0.5 0.5]); 
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    ylabel('Z-Score')
    xlabel('Time (min)')
    setFigureDefaults;
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 6 6])
end

% Plot the autocorrelation of factors 
figure('position',[1105 524 433 454]); hold on;
for i = 1:num_features
    %plot pdf 
    [xc,lags] = xcorr(features_downsampled(:,i)-nanmean(features_downsampled(:,i)),120*13,'coeff'); 
    idx = [ceil(numel(lags)/2)+1:numel(lags)];
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',col(i,:));   
    title({'Behavioral Feature';'Autocorrelation'},'FontName','Arial','FontSize',16,'FontWeight','normal')    
    ylim([min(xc) 1]);
    xlim([0 max(lags/13)])
    ylabel('Rho')
    xlabel('time (s)')
    
    %get halflife
    tau=find(xc(idx)>=0.5*xc(idx(1)),1,'last')/13;
    text(30,0.3+(i*0.1),sprintf('%s \\tau = %.3g s',labels_abbrev{i},tau),'FontSize',16,'FontName','Arial')
    setFigureDefaults;       
end
pos = get(gca,'position');
set(gca,'position',[3 3 6 6])


% Correlation between factors
figure; hold on; 
rho = corr(features_downsampled);
imagesc(rho,[-0.5 1]);
c = colorbar;
ylabel(c,'Correlation','FontSize',16,'FontName','Arial')
colormap magma
set(gca,'XTick',(1:num_features),'YTick',(1:num_features),'YTickLabel',labels_abbrev,'XTickLabel',labels_abbrev,'XTickLabelRotation',45)
xlim([0.5,num_features+0.5])
ylim([0.5,num_features+0.5])
axis square
title({'Correlation Between','Behavioral Features'},'FontSize',16,'FontName','Arial','FontWeight','normal')
setFigureDefaults; 
set(gca,'position',[2 2 6 6])

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','BehavioralCluster_Metrics',fn_path,1);
    close all;
end

%% Tabulate clusters
[clusters,~,indx_clusters] = unique(features_binned,'rows');
temp = [];
unique_states = unique(indx_clusters);
for i = 1:numel(unique_states)
    temp(i) = sum(indx_clusters==unique_states(i));  
end
clusters = clusters(unique_states,:);
clusters(:,end+1) = temp;

%reorder states by decreasing contribution
[~, idx] = sort(temp,'descend');
clusters = clusters(idx,:);

%rename indx_clusters
indx_clusters_temp = NaN(size(indx_clusters));
for i = 1:numel(idx)
    indx_clusters_temp(indx_clusters==idx(i))=i;
end   
indx_clusters = indx_clusters_temp;
unique_states = unique(indx_clusters);

%create a 'junk' state that groups together all states that occur for less 1% of activity
bad_states = unique_states(clusters(:,num_features+1)<=ceil(1/100*numel(indx_clusters)));
if ~isempty(bad_states)
    bad_indices = ismember(indx_clusters,bad_states);

    %combine residual states
    unique_states(ismember(unique_states,bad_states))=bad_states(1);
    indx_clusters(ismember(indx_clusters,bad_states))=bad_states(1);

    % indx_clusters = movmode(indx_clusters,bp.movmode_dur);
    temp = [];
    unique_states = unique(indx_clusters);
    for i = 1:numel(unique_states)
        temp(i) = sum(indx_clusters==unique_states(i));  
    end
    clusters = clusters(unique_states,:);
    clusters(:,end) = temp;
end
num_states = numel(clusters(:,1));


%% shuffle the order of each behavioral state

%break each motif state into a set of consecutive indices
frame_indices = arrayfun(@(x) find(indx_clusters==x),unique_states,'UniformOutput',0); 
consecutive_frames = cellfun(@(x) find(diff([false;[1;diff(x)]==1;false])~=0), frame_indices,'UniformOutput',0);
consecutive_frames = cellfun(@(x) reshape(x', 2,[])', consecutive_frames,'UniformOutput',0);

%break into cell array of cells with the different consecutive indices
cluster_groups = {};
shuffled_labels_organized = {};
for j = 1:numel(frame_indices)
   x = frame_indices{j};
   y = consecutive_frames{j};
   for i = 1:size(y,1)
       if i == 1
           cluster_groups{j,i} = x(1:y(i,2)-1);
       elseif i == size(y,1)
           cluster_groups{j,i} = x(y(i-1,2):end);
       else
           cluster_groups{j,i} = x(y(i-1,2):y(i,2)-1);
       end
       shuffled_labels_organized{j,i} = ones(1,numel(cluster_groups{j,i}))*unique_states(j);
   end
end
shuffled_labels = shuffled_labels_organized(:);

%%
num_features = size(features_downsampled,2);
n_shuf = 1000;
rng('default')

[behav_state_pev, state_avg] = TrialPEV(features_downsampled,indx_clusters);

behav_state_pev_shuf = NaN(n_shuf,size(features_downsampled,2));
for cur_shuf = 1:n_shuf    
   temp = cat(2,shuffled_labels{randperm(numel(shuffled_labels))}); 
   [behav_state_pev_shuf(cur_shuf,:), ~] = TrialPEV(features_downsampled,temp);  
end

% add the original to the shuffle
behav_state_pev_shuf = cat(1,behav_state_pev_shuf,behav_state_pev);

%pval 
pval = sum(behav_state_pev_shuf>=behav_state_pev)/size(behav_state_pev_shuf,1);

% plot permutation test
figure('position',[215 116 1259 862]); hold on; 
[r,c] = numSubplot(size(behav_state_pev,2),1);
for cur_feature = 1:num_features
   subplot(r,c,cur_feature);
   histogram(behav_state_pev_shuf(:,cur_feature)*100,'BinWidth',0.1,'EdgeColor','none','FaceColor',[0.5 0.5 0.5]);
   line([behav_state_pev(:,cur_feature)*100,behav_state_pev(:,cur_feature)*100],get(gca,'ylim'),'linestyle','--','color','r','linewidth',2);
   ylabel('Count');
   xlabel('Percent Explained Variance');
   title(sprintf('State %s Captures Significant\nVariance in Behavioral State',labels_abbrev{cur_feature}),'FontName','Arial','FontSize',16,'FontWeight','normal');
   
   %add pvalue
   text(max(get(gca,'xlim'))/2, max(get(gca,'ylim'))/2, sprintf('p = %.2g',pval(cur_feature)),'FontName','Arial','FontSize',16);
   
   setFigureDefaults;
   pos = get(gca,'position');
   set(gca,'position',[pos(1) pos(2) 6 6]);   
end

% expvar bar plot
figure; hold on; 
bar(behav_state_pev*100,'edgecolor','k','facecolor',[0.75 0.75 0.75],'LineWidth',1);
set(gca,'XTick',(1:numel(pval)),'XTickLabels',labels_abbrev,'XTickLabelRotation',90)
for i = 1:numel(pval)
    AddSig(1,pval(i),[i i behav_state_pev(i)*100 behav_state_pev(i)*100],4,15,1,45)
end
ylabel('Percent Explained Variance');
title({'Behavioral States Capture';'Significant Variance in';'Measured Features'},'FontName','Arial','FontSize',16,'Fontweight','normal')
setFigureDefaults
set(gca,'position',[3 3 6 6])
ylim([0 100])

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Permutation_Test_BehavioralFeatures',fn_path,1);
    close all;
end

%% 
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:); 
H = data.H(expvaridx,:);

H_weight = NaN(size(H))';
for cur_motif = 1:size(H,1)
    H_weight(:,cur_motif) = helper.reconstruct(nanmean(w(:,cur_motif,:),1),H(cur_motif,:)); 
end %motif loop


%% each time the state occurs, get the average intensity of each motif
temp = cellfun(@(x) numel(x)>1, cluster_groups,'UniformOutput',0);
cluster_group_clean = cluster_groups;
for i = 1:size(temp,1)
    cluster_group_clean(i,cell2mat(temp(i,:))==0)={[]};
end

%loop through each cluster and get the average weighting of each motif per 
avg_h = [];
avg_pev = {};
for i = 1:size(cluster_group_clean,1)
   temp = cluster_group_clean(i,:);
   temp = temp(~cellfun('isempty',temp));
   temp = cellfun(@(x) SplitIntoSubarray(x,7),temp,'UniformOutput',0);
   temp = temp(~cellfun('isempty',temp));
   temp = cat(2,temp{:});  
   temp = temp(:,1:2:end);
   for j = 1:size(temp,2)
      temp_log = H_weight(temp(:,j),:);
      temp_log(temp_log<=eps)=NaN;
      avg_h{i}(j,:) = nanmean(log(temp_log)); 
%       avg_h{i}(j,:) = H_weight(temp(:,j),:);
   end
end

if ~isempty(bad_states)
    avg_h(num_states:end)=[];
end

if ~isempty(bad_states)
    state_code=clusters(1:end-1,:);
else
    state_code=clusters;
end

% seperate samples by 0.5 second to get more independent samples
% avg_h = cellfun(@(x) x(1:7:end,:), avg_h, 'UniformOutput',0);

%%
%get distribution
temp_avg = MakeCellsEqual(avg_h,1,1); 
temp_avg = cat(3,temp_avg{:});
temp_avg(temp_avg<=eps)=NaN;
temp_avg(isnan(temp_avg))=[];

figure; hold on; 
[f,xi] = ksdensity((temp_avg(:))); 
plot(xi,f,'linewidth',2,'color',[0.25 0.25 0.25]);
title({'Distribution of Motif';'Intensity Before log Transform'},'FontName','Arial','FontSize',16,'FontWeight','normal');
xlabel('Motif Intensity');
ylabel('PDF');
setFigureDefaults
set(gca,'position',[3 3 6 6]);

figure; hold on; 
[f,xi] = ksdensity(log(temp_avg(:))); 
plot(xi,f,'linewidth',2,'color',[0.25 0.25 0.25]);
title({'Distribution of Motif';'Intensity After log Transform'},'FontName','Arial','FontSize',16,'FontWeight','normal');
xlabel('ln(motif intensity)');
ylabel('PDF');
setFigureDefaults
set(gca,'position',[3 3 6 6]);

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Distribution_figures',fn_path,1);
    close all;
end


%% for each motif, get the average +/- the sem. anova comparing that motif across behavioral states
fp = fig_params;
sig_color = [0.75 0.75 0.75; 0.25 0.25 0.25];
col = getColorPalet(numel(unique_states));
temp_avg = MakeCellsEqual(avg_h,1,1); 
temp_avg = cat(3,temp_avg{:});
pval = [];
fstat = [];
sig_comparison = {};
for i = 1:size(temp_avg,2)
   test = squeeze(temp_avg(:,i,:));
   group = ones(size(test)).*(1:size(test,2));
   test = test(:);
%    test(test<=eps)=NaN;
   group = group(:);
   group(isnan(test))=[];
   test(isnan(test))=[];     
   [pval(i),tab,stats] = anova1(log(test),group,'off');
   c = multcompare(stats,'Display','off');
   temp =c(:,1:2);
   sig_comparison{i} = temp(sum(c(:,end)<0.05,2)>0,:);   
   fstat(i) = tab{2,5};
end

figure('Position',[0 0 1000 1000]); hold on; 
s3 = subplot(313,'Units','centimeters','Position',[4 6 2  4.5]); hold on
s1 = subplot(312,'Units','centimeters','Position',[8 6 10 4.5]); hold on
s2 = subplot(311,'Units','centimeters','Position',[8 11 10 1.5]); hold on

axes(s1);
y=[]; 
hold on
for cur_state = 1:size(temp_avg,3)
    for cur_motif = 1:size(temp_avg,2)
        temp = squeeze(temp_avg(:,cur_motif,cur_state));
%         temp(temp<=eps)=NaN;        
        y(cur_state,cur_motif) = nanmean((temp));
    end    
end
imagesc(y./nanmean(y),[0.95 1.05])
colormap(gca,flipud(redgreencmap(64,'Interpolation','linear')));
c = colorbar;
set(c,'YTick',[0.95, 1, 1.05]);
ylabel(c,{'Normalized Motif';'Intensity (log)'},'FontSize',16,'Fontweight','normal','FontName','Arial');
set(c,'units','centimeters','position',[18.25 6 0.5 4.5])

for x_grid = 0.5:1:size(y,2)+0.5
    line([x_grid,x_grid],[0.5,size(y,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(y,1)+0.5
    line([0.5,size(y,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
xlabel('Basis Motifs')

xlim([0.5 size(y,2)+.5])
ylim([0.5 size(y,1)+0.5])
set(gca,'YColor','none')
setFigureDefaults
% For each significant comparison put a box around motifs that were different
for i = 1:numel(sig_comparison)
   temp = unique(sig_comparison{i}(:));

   for j = 1:numel(temp)
       line([i-0.5, i+0.5],[temp(j)-0.5, temp(j)-0.5],'color',[0.2 0.7 1],'linewidth',2.25)
       line([i-0.5, i+0.5],[temp(j)+0.5, temp(j)+0.5],'color',[0.2 0.7 1],'linewidth',2.25)
       line([i-0.5, i-0.5],[temp(j)-0.5, temp(j)+0.5],'color',[0.2 0.7 1],'linewidth',2.25)
       line([i+0.5, i+0.5],[temp(j)-0.5, temp(j)+0.5],'color',[0.2 0.7 1],'linewidth',2.25)
   end
    
end

axes(s2); hold on
%Plot the fscore and the significance
plot(fstat,'LineWidth',2,'Marker','.','MarkerSize',5,'MarkerEdgeColor',[0.4 0.4 0.4],'Color',[0.4 0.4 0.4])
for i = 1:numel(pval)
   AddSig(1,pval(i),[i-0.1,i-0.1,fstat(i),fstat(i)],1,5,1,90)
end
%Change marker color for significant motifs (lighter)
scatter((1:1:size(temp_avg,2)),fstat,50,sig_color((pval<=0.05/numel(pval))+1,:),'filled')
xlim([0.5 size(temp_avg,2)+.5])
set(gca,'XTick',(1:2:size(temp_avg,2)),'Ytick',[min(get(gca,'Ytick')),max(get(gca,'Ytick'))]);
ylabel({'F';'stat';''},'Rotation',0,'Units','Centimeters','position',[11.25 1.4]);
set(gca,'yaxislocation','right')
set(gca,'XColor','none');
setFigureDefaults

axes(s3); hold on
imagesc(state_code(:,1:num_features),[0 1.25]); 
colormap(gca,'magma')
ylabel('Behavioral State')
set(gca,'XTick',(1:num_features),'XTickLabel',labels,'XTickLabelRotation',90,'TickLength',[0,0],'YTick',(1:1:size(state_code,1)))
set(gca,'Xaxislocation','top')
for x_grid = 0.5:1:num_features+0.5
    line([x_grid,x_grid],[0.5,size(state_code,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(state_code,1)+0.5
    line([0.5,num_features+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
for i = 1:size(state_code,1)
   text(num_features+1,i,sprintf('%.2g%%',sum(indx_clusters==unique_states(i))/numel(indx_clusters)*100),'FontName','Arial','FontSize',16,'FontWeight','normal');
end
xlim([0.5 num_features+0.5])
ylim([0.5 size(state_code,1)+0.5])
setFigureDefaults

set(gcf,'Position',[680   150  875  650]);

fh = gcf;

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','ANOVA_figure',fn_path,1);
    close all;
    save([fn_path, 'indiviudalAnova.mat'],'sig_comparison','fstat','pval');
end
%% Plot the H autocorrelation
figure('position',[ 209         101        1329         877]); hold on; 
H_weight_temp = (H_weight-nanmean(H_weight,1))';
tau = NaN(1,size(H_weight_temp,1));
% [r,c] = numSubplot(size(H_weight_temp,1),1);
for i = 1:size(H_weight_temp,1)
    %plot pdf     
%     subplot(r,c,i); hold on; 
    [xc,lags] = xcorr(H_weight_temp(i,:),10*13,'coeff'); 
    
    N = (numel(H_weight_temp(i,:)));
    for j = 1:numel(xc)
        t = xc(j)*(sqrt(N-2)/sqrt((1-xc(j)^2)));        
%         pval(j) = 1-tcdf(t,N-2); %right tailed stat
        s = tcdf(t,N-2);
        pval(j) = 2 * min(s,1-s);
    end
    
    idx = [ceil(numel(lags)/2)+1:numel(lags)];
    pval = pval(idx);
    sig_xc = xc(idx);
    sig_xc(pval>=0.05/numel(idx))=NaN;    
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',[0.2 0.2 0.2]);
    plot(lags(idx)/13,sig_xc,'linewidth',2,'color',[0.8500 0.3723 0.0078]);
    
    
    title(sprintf('Motif %d',i),'FontName','Arial','FontSize',12,'FontWeight','normal')    
    xlim([0 max(lags/13)])
    ylabel('Rho')
    xlabel('time (s)')    
    tau(i)=find(xc(idx)>=0.5*xc(idx(1)),1,'last')/13;
    
    line([0 max(lags/13)],[0 0],'linestyle','--','linewidth',2,'color',[0.75 0.75 0.75]); 
    text(3,0.5,sprintf('\\tau = %.2gs',tau(i)),'FontSize',12,'FontName','Arial')
    setFigureDefaults; 
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 3 3])
    
    
    ylim([-0.25 1]);
end

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','H_autocorrelation_sigtesting',fn_path,1);
    close all;
end


%% save off consecutive example frames (max 5 instances) for each cluster
frame_indices = arrayfun(@(x) find(indx_clusters==x),unique_states,'UniformOutput',0); 
% get starts and stops of consecutive indices
consecutive_frames = cellfun(@(x) find(diff([false;[1;diff(x)]==1;false])~=0), frame_indices,'UniformOutput',0);
consecutive_frames = cellfun(@(x) reshape(x', 2,[])', consecutive_frames,'UniformOutput',0);
% select up to 5 longest instances
[t, temp] = cellfun(@(x) maxk(x(:,2)-x(:,1),min([5,size(x,1)])),consecutive_frames,'UniformOutput',0);
consecutive_frames =  cellfun(@(x,y) x(y,:),  consecutive_frames, temp, 'UniformOutput',0);

% Plot ethogram
figure('position',[462 549 1029 351]); hold on; 
x = linspace(0,60,numel(indx_clusters));
for i = 1:numel(frame_indices)
   plot(x(frame_indices{i}),ones(1,numel(frame_indices{i}))*i,'linestyle','none','marker','.','markersize',5,'color','k')
end
ylim([0 numel(frame_indices)+1]);
xlabel('Time (min)')
ylabel('Behavioral State')
set(gca,'TickLength',[0,0])
title('Ethogram of Behavioral States','FontName','Arial','Fontsize',16,'FontWeight','normal')
setFigureDefaults
set(gca,'position',[3 3 18.5 5])    
%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-png','Ethogram',fn_path,1);
    close all;
end

%% Now get the movie frames in the origin film that correspond to each cluster and save off those snippets
bodycam = VideoReader([fn_path fn_bodycam]);
facecam = VideoReader([fn_path fn_facecam]);

for i = 1:numel(consecutive_frames)
    savedir = [fn_path sprintf('BehavioralCluster_%d',i)];
    if ~exist(savedir)
        mkdir(savedir)
    end

    %take shorter if the clip is super long
    for cur_snip = 1:size(consecutive_frames{i},1)
       %get the frame indices in the original videos
       temp = frame_indices{i}(consecutive_frames{i}(cur_snip,1):(consecutive_frames{i}(cur_snip,2)-1));    
       
       %get the coordinates in the original space + the ONSET
       temp = [min(temp), max(temp)];
       temp = round(x_query_ds(temp),0);
       temp = temp+onset;
       
       if numel(temp)~=2
           continue
       end
       
       
       if temp(2)-temp(1)<60
           continue
       end 
       
       if temp(2)-temp(1)>(15*60)
          temp(2)=temp(1)+(15*60);
       end
                    
       
       %add 1 second on either side
       temp(1) = temp(1)-(60);
       temp(2) = temp(2)+(60);
       
       v = VideoWriter([savedir filesep sprintf('BodyCamSnippet%d',cur_snip)],'Motion JPEG AVI');
       v.FrameRate = 60; 
       open(v)
       frames = read(bodycam,temp);
       writeVideo(v,frames);   
       close(v)
       
       v = VideoWriter([savedir filesep sprintf('FaceCamSnippet%d',cur_snip)],'Motion JPEG AVI');
       v.FrameRate = 60; 
       open(v)
       frames = read(facecam,temp);
       writeVideo(v,frames);   
       close(v)
       
    end   
end

%% Create inidviudal high res frames for select snippets

cur_cluster = 3;
cur_snip = 2; 
savedir = [fn_path sprintf('BehavioralCluster_%d_Frames',cur_cluster)];
if ~exist(savedir)
    mkdir(savedir)
end

temp = frame_indices{cur_cluster}(consecutive_frames{cur_cluster}(cur_snip,1):(consecutive_frames{cur_cluster}(cur_snip,2)-1));    

%get the coordinates in the original space + the ONSET
temp = [min(temp), max(temp)];
temp = round(x_query_ds(temp),0);
temp = temp+onset;

if temp(2)-temp(1)>(15*60)
  temp(2)=temp(1)+(15*60);
end

frames = read(bodycam,temp);
%just sample 4x a second
for cur_frame = 1:15:size(frames,4)
   figure; 
   imshow(frames(100:325,175:500,:,cur_frame))
   axis off           
   setFigureDefaults
   set(gca,'position',[3 3 6 6])
   drawnow
end
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg',sprintf('BodySnippet%d',cur_snip),savedir,1);
close all

frames = read(facecam,temp);
%just sample 4x a second
for cur_frame = 1:15:size(frames,4)
   figure; 
   imshow(frames(:,:,:,cur_frame))
   axis off           
   setFigureDefaults
   set(gca,'position',[3 3 6 6])   
   drawnow
   
end
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg',sprintf('FaceCamSnippet%d',cur_snip),savedir,1);
close all

% % % % 
% % % % %%temp_avg = MakeCellsEqual(avg_h,1,1); 
% % % % temp_avg = cat(3,temp_avg{:});
% % % % temp_avg(temp_avg<=eps)
% % % % motif_group = ones(size(temp_avg)).*(1:size(temp_avg,2));
% % % % state_group = temp_avg; 
% % % % for i = 1:size(temp_avg,3)
% % % %    state_group(:,:,i) = ones(size(state_group(:,:,i)))*i; 
% % % % end
% % % % temp_avg(temp_avg<=eps)=NaN;
% % % % temp_avg = log(temp_avg);
% % % % temp_avg = temp_avg(:);
% % % % motif_group = motif_group(:);
% % % % state_group = state_group(:);
% % % % [p, tbl,stats] = anovan(temp_avg(1:7:end),cat(2,motif_group(1:7:end),state_group(1:7:end)),'model','interaction');
% % % % %%
% % % % save([fn_path '2wayAnovaData.mat'],'p','tbl','stats');