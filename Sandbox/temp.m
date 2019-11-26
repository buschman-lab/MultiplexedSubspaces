%add paths
% addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\'));
%set filepaths
fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse431_10_17_2019\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_bodycam = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
fn_dlc = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse431_10_17_2019';
fn_widefield = '431-10-17-2019_1Fitted_block_hemoflag0_1';

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

H = load([fn_path, fn_widefield],'H');
num_frames = size(H.H,2);

savefigs = 1; 

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

%get the distance between nose and front paws
% [frontpaws_to_nose, ~, ~] = parse_dlc(raw_data,{'frontrightpawcenter','frontleftpawcenter'},'nosetip',bp.dlc_epsilon);
% %you've subtracted out the reference. Now get euclidean distance
% frontpaws_to_nose = sqrt(sum(frontpaws_to_nose.^2,2));

%get the distance from nose to tail
% [tail_to_nose, ~, ~] = parse_dlc(raw_data,{'tailroot'},'nosetip',bp.dlc_epsilon);
% tail_to_nose = sqrt(sum(tail_to_nose.^2,2));

%combined all features
% face_motion_energy{2} = -1 * face_motion_energy{2}+max(face_motion_energy{2}(:)); %may need to flip the whisker energy if high whisking actually blurs the camera and make low energy  
features = cat(2,face_motion_energy{:},limb_speed);
% labels = {'nose motion energy','whisker motion energy','front paws to nose distance','limb speed'};
% labels_abbrev = {'NME','WME','F2N','LS'};

labels = {'nose motion energy','whisker motion energy','limb speed'};
labels_abbrev = {'NME','WME','LS'};

features_no_smooth = features;
for i = 1:size(features,2)
    features(:,i) = convn(features(:,i),ones(65,1)/65,'same');
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


%% Phenograph
[ovr_comm_resampled, ovr_Q] = PhenographSimple(features_downsampled(1:3:end,:), 'knn', 15);

%%

%% Plot example behavioral traces
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
%% Plot statistics about the behavioral traces
% split_idx = {0,1,-2.125,.25}; %mouse 432 10_17_2019
% split_idx = {-0.315,0.907,-2.484,-.14}; %mouse 431;
% split_idx = {-0.834,.98,-1,.52}; %mouse 494 10_17_2019

split_idx = {11.96,10.98,0.29};

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
    text(30,0.3+(i*0.1),sprintf('%s \\tau = %.2g s',labels_abbrev{i},tau),'FontSize',16,'FontName','Arial')
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

%% Parse states
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

%% Filter and Tabulate clusters
[clusters,~,indx_clusters] = unique(features_binned,'rows');

% indx_clusters = movmode(indx_clusters,bp.movmode_dur);
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

%remove any states that occur for less 1% of activity
bad_states = unique_states(clusters(:,4)<=ceil(1/100*numel(indx_clusters)));
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
clusters(:,end+1) = temp;


%% shuffle the order of each behavioral state

%break each motif state into a set of consecutive indices
frame_indices = arrayfun(@(x) find(indx_clusters==x),unique_states,'UniformOutput',0); 
consecutive_frames = cellfun(@(x) find(diff([false;[1;diff(x)]==1;false])~=0), frame_indices,'UniformOutput',0);
consecutive_frames = cellfun(@(x) reshape(x', 2,[])', consecutive_frames,'UniformOutput',0);

%break into cell array of cells with the different consecutive indices
cluster_groups = {};
shuffled_labels = {};
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
       shuffled_labels{j,i} = ones(1,numel(cluster_groups{j,i}))*unique_states(j);
   end
end
shuffled_labels = shuffled_labels(:);

%%
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
n_shuf=1000;
rng('default')
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:); 
H = data.H(expvaridx,:);

H_weight = NaN(size(H))';
for cur_motif = 1:size(H,1)
    H_weight(:,cur_motif) = helper.reconstruct(nanmean(w(:,cur_motif,:),1),H(cur_motif,:)); 
end %motif loop

temp = H_weight./squeeze(nanmean(w,[1,3]));
[H_thresh, ~, ~] = ThresholdMatrix(temp,0.15);
H_thresh(H_thresh==0)=NaN;
H_thresh = H_thresh(1:3:end,:);

fprintf('Fraction time active per motif')
fprintf('\n%.2g', sum(H_thresh>0)/size(H_thresh,1));

fprintf('Total time active across motifs: %0.2g',sum(sum(H_thresh>0,2)>0)/size(H_thresh,1))

[motif_pev, ~] = TrialPEV(H_thresh,indx_clusters);

motif_pev_shuf = NaN(n_shuf,size(H_thresh,2));
% motif_avg_shuf = NaN(size(motif_avg,1),size(motif_avg,2),n_shuf);
for cur_shuf = 1:n_shuf    
   temp = cat(2,shuffled_labels{randperm(numel(shuffled_labels))}); 
%    temp = indx_clusters(randperm(numel(indx_clusters)));
   [motif_pev_shuf(cur_shuf,:), ~] = TrialPEV(H_thresh,temp);  
end

% add the original to the shuffle
motif_pev_shuf = cat(1,motif_pev_shuf,motif_pev);
% motif_avg_shuf = cat(3,motif_avg_shuf,motif_avg);

%pval 
pval = sum(motif_pev_shuf>motif_pev)/size(motif_pev_shuf,1);
% pval_avg = sum(motif_avg_shuf>=repmat(motif_avg,1,1,size(motif_avg_shuf,3)),3)/size(motif_avg_shuf,3);
% pval_bin = 1-binocdf(sum(pval_avg<=0.1),size(pval_avg,2),0.1);

% plot pev bar
figure; hold on; 
bar(motif_pev*100,'edgecolor','k','facecolor',[0.75 0.75 0.75],'LineWidth',1);
set(gca,'XTick',(1:2:numel(pval)))
xlabel('Basis Motif');
for i = 1:numel(pval)
    AddSig(1,pval(i),[i i motif_pev(i)*100 motif_pev(i)*100],4,2.5,1,90)
end
ylabel('Percent Explained Variance');
title({'Behavioral States Capture';'Significant Variance in';'Measured Features'},'FontName','Arial','FontSize',16,'Fontweight','normal')
setFigureDefaults
set(gca,'position',[3 3 6 6])
ylim([0 30])


%% Get the percent active during each cluster
temp = [];
for i = 1:numel(unique_states)
    for j = 1:size(H_thresh,2)
        temp(i,j) = sum(H_thresh(indx_clusters==unique_states(i),j)>0)/numel(H_thresh(indx_clusters==unique_states(i),j))*100;
    end
end

figure('position',[680   420   560   558]); hold on; 
imagesc(temp',[0 15])
colormap(magma)
c = colorbar;
ylabel(c,{'% Time Active'},'FontSize',16)
ylim([0.5,size(motif_avg,2)+0.5]);
xlim([0.5,size(motif_avg,1)+0.5]);
xlabel('Behavioral State')
ylabel('Basis Motif')
setFigureDefaults;
set(gca,'position',[3 4 8 8],'box','on');
set(c,'units','centimeters','position',[11.5 4 0.5 8]);
for x_grid = 0.5:1:size(motif_avg,1)+0.5
line([x_grid,x_grid],[0.5,size(motif_avg,2)+0.5],'linewidth',1.5,'color','k')
end
for y_grid = 0.5:1:size(motif_avg,2)+0.5
line([0.5,size(motif_avg,1)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')
end


%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Permutation_Test_H_weightings',fn_path,1);
    close all;
end

%% Plot H distribution 
col = getColorPalet(14);
figure('position',[680   101   858   877]); hold on; 
for i = 1:size(H_weight,2)
    %plot pdf 
%     temp_w = squeeze(nanmean(w,[1,3]));
%     temp = H_weight(:,i)./temp_w(i);
    temp= H_thresh(:,i);
    [f,xi] = ksdensity(temp); 
    plot(xi,f,'linewidth',2,'color',col(i,:)); 
    title(sprintf('Motif %d',i),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;
    ylabel('PDF')
    xlabel('Intensity')
    yvals = get(gca,'ylim');    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) 3 3])
end

%% Plot the Xhat distribution
temp = load(fn_widefield,'data_test');
figure('position',[680   101   858   877]); hold on; 
%plot pdf 
temp = nanmean(temp.data_test,1);
[f,xi] = ksdensity(temp); 
plot(xi,f,'linewidth',2,'color',[0.5 0.5 0.5]); 
title(sprintf('Distribution of average pixel intensity across active pixels %d',i),'FontName','Arial','FontSize',16,'FontWeight','normal')
setFigureDefaults;
ylabel('PDF')
xlabel('Intensity')
yvals = get(gca,'ylim');    
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 3 3])

fprintf('Percent Frames Active = %0.2g', sum(temp>=0.001)/numel(temp)*100)

%% Plot the H autocorrelation
figure('position',[680   101   858   877]); hold on; 
H_weight_temp = (H_weight-nanmean(H_weight,1))';
tau = NaN(1,size(H_weight_temp,1));
for i = 1:size(H_weight_temp,1)
    %plot pdf     
    [xc,lags] = xcorr(H_weight_temp(i,:),20*13,'coeff'); 
    idx = [ceil(numel(lags)/2)+1:numel(lags)];
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',[0.5 0.5 0.5]);   
    title({'Autocorrelation of Motif';'Temporal Weightings'},'FontName','Arial','FontSize',16,'FontWeight','normal')    
    xlim([0 max(lags/13)])
    ylabel('Rho')
    xlabel('time (s)')    
    tau(i)=find(xc(idx)>=0.5*xc(idx(1)),1,'last')/13;
end

% ylim([-0.5 1]);
%get halflife
text(3,0.5,sprintf('\\tau = %.2g +/- %.2gs',nanmean(tau),sem(tau,2)),'FontSize',16,'FontName','Arial')
setFigureDefaults; 

pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 6 6])
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','H_autocorrelation',fn_path,1);
    close all;
end

%% Plot the binary state ID in vertical format
close all; figure('position',[100 50 700, 900]); 
imagesc(clusters(:,1:(end-1)),[0 1.25]); 
colormap magma
ylabel('Behavioral State')

set(gca,'XTick',(1:numel(labels)),'XTickLabel',labels,'XTickLabelRotation',90,'TickLength',[0,0],'YTick',(1:1:size(clusters,1)))
for x_grid = 0.5:1:numel(labels)+0.5
    line([x_grid,x_grid],[0.5,size(clusters,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(clusters,1)+0.5
    line([0.5,numel(labels)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
for i = 1:numel(unique_states)
   text(5,i,sprintf('%.2g%%',sum(indx_clusters==unique_states(i))/numel(indx_clusters)*100),'FontName','Arial','FontSize',16,'FontWeight','normal');
end
text(5.5,15.5,{'Contribution'},'horizontalalignment','center','FontName','Arial','FontSize',16,'FontWeight','normal','Rotation',90)
setFigureDefaults
set(gca,'position',[4 8 3.75 12],'box','on')


if savefigs
    handles = get(groot, 'Children');
    saveas(handles,[fn_path 'BinaryState_Vertical'],'svg')
    close all;
end



%% save off consecutive example frames (allow up to 5 missed frames) (max 5 instances) for each cluster
frame_indices = arrayfun(@(x) find(indx_clusters==x),unique_states,'UniformOutput',0); 
% get starts and stops of consecutive indices
consecutive_frames = cellfun(@(x) find(diff([false;[1;diff(x)]==1;false])~=0), frame_indices,'UniformOutput',0);
consecutive_frames = cellfun(@(x) reshape(x', 2,[])', consecutive_frames,'UniformOutput',0);
% select up to 5 longest instances
[t, temp] = cellfun(@(x) maxk(x(:,2)-x(:,1),min([5,size(x,1)])),consecutive_frames,'UniformOutput',0);
consecutive_frames =  cellfun(@(x,y) x(y,:),  consecutive_frames, temp, 'UniformOutput',0);

%% Plot ethogram
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
%%
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


