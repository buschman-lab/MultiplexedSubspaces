%add paths
% addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\'));
%set filepaths
fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse494_10_17_2019\';
fn_facecam = 'Cam_0_20191017-184642.avi';
fn_bodycam = 'Cam_1_20191017-184642_Mouse494_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
fn_dlc = 'Cam_1_20191017-184642_Mouse494_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse494_10_17_2019';
fn_widefield = '494-10-17-2019_1Fitted_block_hemoflag0_1';

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
[frontpaws_to_nose, ~, ~] = parse_dlc(raw_data,{'frontrightpawcenter','frontleftpawcenter'},'nosetip',bp.dlc_epsilon);
%you've subtracted out the reference. Now get euclidean distance
frontpaws_to_nose = sqrt(sum(frontpaws_to_nose.^2,2));

%get the distance from nose to tail
[tail_to_nose, ~, ~] = parse_dlc(raw_data,{'tailroot'},'nosetip',bp.dlc_epsilon);
tail_to_nose = sqrt(sum(tail_to_nose.^2,2));

%combined all features
face_motion_energy{2} = -1 * face_motion_energy{2}; %may need to flip the whisker energy if high whisking actually blurs the camera and make low energy  
features = cat(2,face_motion_energy{:},frontpaws_to_nose,limb_speed);
labels = {'nose motion energy','whisker motion energy','front paws to nose distance','limb speed'};
labels_abbrev = {'NME','WME','F2N','LS'};

features_no_smooth = features;
for i = 1:size(features,2)
    features(:,i) = convn(features(:,i),ones(65,1)/65,'same');
end

%Trim to match start and stop of imaging 
features = features(onset:offset,:);
features_no_smooth = features_no_smooth(onset:offset,:);

% Downsample to match motif duration
features_downsampled = NaN(num_frames,size(features,2));
features_no_smooth_downsampled = NaN(num_frames,size(features_no_smooth,2));
for i = 1:size(features_downsampled,2)
    temp = features(:,i);
    features_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');

    temp = features_no_smooth(:,i);
    features_no_smooth_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');    
end

%store the mapping from original to downsampled 
x_query_ds = linspace(1,size(features,1),num_frames);

if bp.zscore %optional zscore
   features_downsampled = zscore(features_downsampled,1);
   features_no_smooth_downsampled = zscore(features_no_smooth_downsampled,1);
end

clear raw_data facecam_data W_clust_smooth

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
split_idx = {-0.834,.98,-1,.52}; %mouse 494 10_17_2019

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
    [xc,lags] = xcorr(features_downsampled(:,i),120*13,'coeff'); 
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

% Plot the autocorrelation of the unsmooth factors to confirm appropriate smoothing
figure('position',[1105 524 433 454]); hold on; 
for i = 1:num_features
    %plot pdf 
    [xc,lags] = xcorr(features_no_smooth_downsampled(:,i),120*13,'coeff'); 
    idx = [ceil(numel(lags)/2)+1:numel(lags)];
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',col(i,:));   
    title({'Behavioral Feature';'Autocorrelation (NOSMOOTH)'},'FontName','Arial','FontSize',16,'FontWeight','normal')    
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

%remove any states that occur for less 0.1% of activity
bad_states = unique_states(clusters(:,5)<=ceil(0.1/100*numel(indx_clusters)));
bad_indices = ismember(indx_clusters,bad_states);
clusters(clusters(:,5)<=ceil(0.1/100*numel(indx_clusters)),:)=[];
unique_states(ismember(unique_states,bad_states))=[];
features_downsampled(ismember(indx_clusters,bad_states),:)=[];
features_binned(ismember(indx_clusters,bad_states),:)=[];
indx_clusters(ismember(indx_clusters,bad_states))=[];



warning('you are removing %d states that contrinubte <0.1% of activity',numel(bad_states));

%% Look at the relative frame-wise explained variance of each motif for each state
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:);
H = data.H(expvaridx,:);
H(:,bad_indices)=[];
data.data_test(:,bad_indices)=[];

load_frame = NaN(size(H,1),size(H,2));
for cur_motif = 1:size(H,1)
    fprintf('\nworking on motif %d',cur_motif);
    WH = helper.reconstruct(w(:,cur_motif,:),H(cur_motif,:));  
    for cur_f = 1:size(H,2)
        load_frame(cur_motif,cur_f) = 1 - nanvar(data.data_test(:,cur_f)-WH(:,cur_f))./nanvar(data.data_test(:,cur_f));        
    end   
end %motif loop
clear data WH

%remove motifs that capture less than 1 percent of the variance to a given frame
load_frame(load_frame<0.01)=NaN;

%split by behavioral state
expvar_state = arrayfun(@(x) load_frame(:,indx_clusters==x),unique_states,'UniformOutput',0); 

expvar_state = cellfun(@(x) nanmean(x./nansum(x,1),2), expvar_state, 'UniformOutput',0);
expvar_state = [expvar_state{:}]; 


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
rng('default')

% basal_variance = var(features_downsampled,[],1);
%get the variance of behavioral factor during a state
features_variance = NaN(numel(unique_states),size(features_downsampled,2));
for i = 1:numel(unique_states)
    features_variance(i,:) = var(features_downsampled(indx_clusters==unique_states(i),:),[],1);
end

features_variance_shuf = NaN(numel(unique_states),size(features_downsampled,2),1000);
for num_shuf = 1:5000    
   temp = cat(2,shuffled_labels{randperm(numel(shuffled_labels))});   
   for i = 1:numel(unique_states)
      features_variance_shuf(i,:,num_shuf) = var(features_downsampled(temp==unique_states(i),:),[],1);
   end      
end

%Add correct to the shuffled data set
features_variance_shuf = cat(3,features_variance_shuf,features_variance);
pval_store=NaN(num_features,numel(unique_states));
pval_bin = NaN(num_features,1);
for cur_feature = 1:num_features
    figure('position',[215 116 1259 862]); hold on; 
    [r,c] = numSubplot(size(features_variance,1),1);
    for cur_state = 1:size(features_variance,1)
       subplot(r,c,cur_state); 
       histogram(features_variance_shuf(cur_state,cur_feature,:),'numBins',100);
       line([features_variance(cur_state,cur_feature),features_variance(cur_state,cur_feature)],...
           get(gca,'ylim'),'linestyle','--','color','r','linewidth',2);
       
       %get one tailed pvalues
       if features_variance(cur_state,cur_feature)<=nanmean(features_variance_shuf(cur_state,cur_feature,:))
          pval = sum(features_variance_shuf(cur_state,cur_feature,:)<=features_variance(cur_state,cur_feature));          
       else
          pval = sum(features_variance_shuf(cur_state,cur_feature,:)>=features_variance(cur_state,cur_feature));
       end
       pval = pval/numel(features_variance_shuf(cur_state,cur_feature,:));
       pval_store(cur_feature,cur_state) = pval;
       title(sprintf('%s State %d. \nPval %.2g',labels_abbrev{cur_feature},cur_state,pval),'FontSize',16,'FontWeight','normal');
       xlim([get(gca,'xlim') + [-1,1]])
       setFigureDefaults;
       pos = get(gca,'position');
       set(gca,'position',[pos(1) pos(2) 2.5 2.5])
       
    end
    pval_bin(cur_feature) = 1 - binocdf(sum(pval_store(cur_feature,:)<0.05),numel(pval_store(cur_feature,:)),0.05);
end

%normalize feature variance
features_variance_norm = features_variance;
for i = 1:size(features_variance,2)
    features_variance_norm(:,i) = (features_variance_norm(:,i)-min(features_variance_norm(:,i)))/(max(features_variance_norm(:,i))-min(features_variance_norm(:,i)));
end

figure('position',[680   420   560   558]); hold on; 
imagesc(floor(log10(pval_store)),[-4 -1]);
cmap = magma(4);
cmap = (cmap(1:4,:));
colormap(cmap);
c=colorbar;
ylabel(c,'log_{10}(\itp-value)','FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Behavioral Features Are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
set(gca,'YTick',(1:size(pval_store,1)),'YTickLabels',labels_abbrev)
for i = 1:numel(pval_bin)
   text(numel(unique_states)+1,i,sprintf('%.2g',pval_bin(i)),'FontSize',16,'FontName','Arial','FontWeight','normal')
end
text(numel(unique_states)+2.5, -3,{'Binomial';'Probability'},'Rotation',90,'HorizontalAlignment','Center',...
    'FontSize',16,'FontWeight','normal','FontName','Arial')
setFigureDefaults;
set(gca,'position',[3 4 8 2.5],'box','on');
set(c,'units','centimeters','position',[13.5 4 0.5 2.5]);
for x_grid = 0.5:1:size(features_variance_norm,1)+0.5
    line([x_grid,x_grid],[0.5,size(features_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(features_variance_norm,2)+0.5
    line([0.5,size(features_variance_norm,1)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end 
set(c,'YTick',linspace(-4,-1,5),'YTickLabel',{'','-4','-3','-2','-1'})

% plot the relative increase or decrease of the motif's activity
figure('position',[680   420   560   558]); hold on; 
imagesc(features_variance_norm',[0 1]);
colormap(magma)
c=colorbar;
ylabel(c,{'Variance';'(normalized)'},'FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Behavioral Features Are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
set(gca,'YTick',(1:size(pval_store,1)),'YTickLabels',labels_abbrev)
setFigureDefaults;
set(gca,'position',[3 4 8 2.5],'box','on');
set(c,'units','centimeters','position',[12 4 0.5 2.5]);
for x_grid = 0.5:1:size(features_variance_norm,1)+0.5
    line([x_grid,x_grid],[0.5,size(features_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(features_variance_norm,2)+0.5
    line([0.5,size(features_variance_norm,1)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Permutation_Test_BehavioralFeatures',fn_path,1);
    close all;
end


%% Do the same thing with your H weightings
rng('default')
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:);
H = data.H(expvaridx,:);
H(:,bad_indices)=[];

H_weight = NaN(size(H))';
for cur_motif = 1:size(H,1)
    H_weight(:,cur_motif) = helper.reconstruct(nanmean(w(:,cur_motif,:),1),H(cur_motif,:)); 
end %motif loop
clear data H w

%get the variance of behavioral factor during a state
H_variance = NaN(numel(unique_states),size(H_weight,2));
for i = 1:numel(unique_states)
    H_variance(i,:) = var(H_weight(indx_clusters==unique_states(i),:),[],1);
end

H_variance_shuf = NaN(numel(unique_states),size(H_weight,2),1000);
for num_shuf = 1:5000    
   temp = cat(2,shuffled_labels{randperm(numel(shuffled_labels))});   
   for i = 1:numel(unique_states)
      H_variance_shuf(i,:,num_shuf) = var(H_weight(temp==unique_states(i),:),[],1);
   end      
end
%Add correct to the shuffled data set
H_variance_shuf = cat(3,H_variance_shuf,H_variance);
pval_store=NaN(size(H_variance,2),size(H_variance,1));
pval_bin = NaN(size(H_variance,2),1);
for cur_motif = 1:size(H_variance,2)
    figure('position',[215 116 1259 862]); hold on; 
    [r,c] = numSubplot(size(H_variance,1),1);
    for cur_state = 1:size(H_variance,1)
       subplot(r,c,cur_state); 
       histogram(H_variance_shuf(cur_state,cur_motif,:),'numBins',100);
       line([H_variance(cur_state,cur_motif),H_variance(cur_state,cur_motif)],...
           get(gca,'ylim'),'linestyle','--','color','r','linewidth',2);
       
       %get one tailed pvalues
       if H_variance(cur_state,cur_motif)<=nanmean(H_variance_shuf(cur_state,cur_motif,:))
          pval = sum(H_variance_shuf(cur_state,cur_motif,:)<=H_variance(cur_state,cur_motif));          
       else
          pval = sum(H_variance_shuf(cur_state,cur_motif,:)>=H_variance(cur_state,cur_motif));
       end
       pval = pval/numel(H_variance_shuf(cur_state,cur_motif,:));
       pval_store(cur_motif,cur_state) = pval;
       title(sprintf('Motif %d State %d. \nPval %.2g',cur_motif,cur_state,pval),'FontSize',16,'FontWeight','normal');
%        xlim([get(gca,'xlim') *1.5])
       setFigureDefaults;
       pos = get(gca,'position');
       set(gca,'position',[pos(1) pos(2) 2.5 2.5])
    end   
    
    pval_bin(cur_motif) = 1 - binocdf(sum(pval_store(cur_motif,:)<0.1),numel(pval_store(cur_motif,:)),0.1);
end

%normalize motif variance
H_variance_norm = H_variance;
for i = 1:size(H_variance,2)
    H_variance_norm(:,i) = (H_variance_norm(:,i)-min(H_variance_norm(:,i)))/(max(H_variance_norm(:,i))-min(H_variance_norm(:,i)));
end


figure('position',[680   420   560   558]); hold on; 
imagesc(floor(log10(pval_store)),[-4 -1]);
cmap = magma(4);
cmap = (cmap(1:4,:));
% cmap(1,:) = [1,1,1];
colormap(cmap);
c=colorbar;
ylabel(c,'log_{10}(\itp-value)','FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Basis Motifs are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
ylabel('Basis Motif')
for i = 1:numel(pval_bin)
   text(numel(unique_states)+1,i,sprintf('%.2g',pval_bin(i)),'FontSize',16,'FontName','Arial','FontWeight','normal')
end
text(numel(unique_states)+2.5, -3,{'Binomial';'Probability'},'Rotation',90,'HorizontalAlignment','Center',...
    'FontSize',16,'FontWeight','normal','FontName','Arial')
setFigureDefaults;
set(gca,'position',[3 4 8 8],'box','on');
set(c,'units','centimeters','position',[13.5 4 0.5 8]);
for x_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([x_grid,x_grid],[0.5,size(H_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([0.5,size(H_variance_norm,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end 
set(c,'YTick',linspace(-4,-1,5),'YTickLabel',{'','-4','-3','-2','-1'})

% plot the relative increase or decrease of the motif's activity
figure('position',[680   420   560   558]); hold on; 
imagesc(H_variance_norm',[0 1]);
colormap magma
c=colorbar;
ylabel(c,{'Variance';'(Normalized By Motif)'},'FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Motif are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
ylabel('Basis Motif')
setFigureDefaults;
set(gca,'position',[3 4 8 8],'box','on');
set(c,'units','centimeters','position',[11.5 4 0.5 8]);
for x_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([x_grid,x_grid],[0.5,size(H_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([0.5,size(H_variance_norm,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end 
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Permutation_Test_H_variance',fn_path,1);
    close all;
end

%% Do the same thing with your H occurnaces
[~, h_bin, ~] = ThresholdMatrix(H_weight,1);
basal_rate = sum(h_bin,1)./size(h_bin,1);

%get the variance of behavioral factor during a state
H_variance = NaN(numel(unique_states),size(H_weight,2));
for i = 1:numel(unique_states)
    H_variance(i,:) = sum(h_bin(indx_clusters==unique_states(i),:),1)/sum(indx_clusters==unique_states(i));
    H_variance(i,:) = H_variance(i,:)./basal_rate;
end

H_variance_shuf = NaN(numel(unique_states),size(H_weight,2),1000);
for num_shuf = 1:5000    
   temp = cat(2,shuffled_labels{randperm(numel(shuffled_labels))});   
   for i = 1:numel(unique_states)
        H_variance_shuf(i,:,num_shuf) = sum(h_bin(temp==unique_states(i),:),1)/sum(temp==unique_states(i));
        H_variance_shuf(i,:,num_shuf) = H_variance_shuf(i,:,num_shuf)./basal_rate;
   end      
end
%Add correct to the shuffled data set
H_variance_shuf = cat(3,H_variance_shuf,H_variance);
pval_store=NaN(size(H_variance,2),size(H_variance,1));
pval_bin = NaN(size(H_variance,2),1);
for cur_motif = 1:size(H_variance,2)
    figure('position',[215 116 1259 862]); hold on; 
    [r,c] = numSubplot(size(H_variance,1),1);
    for cur_state = 1:size(H_variance,1)
       subplot(r,c,cur_state); 
       histogram(H_variance_shuf(cur_state,cur_motif,:),'numBins',100);
       line([H_variance(cur_state,cur_motif),H_variance(cur_state,cur_motif)],...
           get(gca,'ylim'),'linestyle','--','color','r','linewidth',2);
       
       %get one tailed pvalues
       if H_variance(cur_state,cur_motif)<=nanmean(H_variance_shuf(cur_state,cur_motif,:))
          pval = sum(H_variance_shuf(cur_state,cur_motif,:)<=H_variance(cur_state,cur_motif));          
       else
          pval = sum(H_variance_shuf(cur_state,cur_motif,:)>=H_variance(cur_state,cur_motif));
       end
       pval = pval/numel(H_variance_shuf(cur_state,cur_motif,:));
       pval_store(cur_motif,cur_state) = pval;
       title(sprintf('Motif %d State %d. \nPval %.2g',cur_motif,cur_state,pval),'FontSize',16,'FontWeight','normal');
%        xlim([get(gca,'xlim') *1.5])
       setFigureDefaults;
       pos = get(gca,'position');
       set(gca,'position',[pos(1) pos(2) 2.5 2.5])
    end   
    
    pval_bin(cur_motif) = 1 - binocdf(sum(pval_store(cur_motif,:)<0.1),numel(pval_store(cur_motif,:)),0.1);

end

%normalize motif variance
H_variance_norm = H_variance;
for i = 1:size(H_variance,2)
    H_variance_norm(:,i) = (H_variance_norm(:,i)-min(H_variance_norm(:,i)))/(max(H_variance_norm(:,i))-min(H_variance_norm(:,i)));
end

figure('position',[680   420   560   558]); hold on; 
imagesc(floor(log10(pval_store)),[-4 -1]);
cmap = magma(4);
cmap = (cmap(1:4,:));
% cmap(1,:) = [1,1,1];
colormap(cmap);
c=colorbar;
ylabel(c,'log_{10}(\itp-value)','FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Basis Motifs are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
ylabel('Basis Motif')
for i = 1:numel(pval_bin)
   text(numel(unique_states)+1,i,sprintf('%.2g',pval_bin(i)),'FontSize',16,'FontName','Arial','FontWeight','normal')
end
text(numel(unique_states)+2.5, -3,{'Binomial';'Probability'},'Rotation',90,'HorizontalAlignment','Center',...
    'FontSize',16,'FontWeight','normal','FontName','Arial')
setFigureDefaults;
set(gca,'position',[3 4 8 8],'box','on');
set(c,'units','centimeters','position',[13.5 4 0.5 8]);
for x_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([x_grid,x_grid],[0.5,size(H_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([0.5,size(H_variance_norm,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end 
set(c,'YTick',linspace(-4,-1,5),'YTickLabel',{'','-4','-3','-2','-1'})

% plot the relative increase or decrease of the motif's activity
figure('position',[680   420   560   558]); hold on; 
imagesc(H_variance_norm',[0 1]);
colormap magma
c=colorbar;
ylabel(c,{'Variance';'(Normalized By Motif)'},'FontName','Arial','FontWeight','normal','Fontsize',16,'Interpreter','tex');
title({'Motif are Differentially';'Engaged Across Behavioral States'},'FontSize',16,'Fontweight','normal','FontName','Arial')
ylim([0.5,size(pval_store,1)+0.5]);
xlim([0.5,size(pval_store,2)+0.5]);
xlabel('Behavioral State')
ylabel('Basis Motif')
setFigureDefaults;
set(gca,'position',[3 4 8 8],'box','on');
set(c,'units','centimeters','position',[11.5 4 0.5 8]);
for x_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([x_grid,x_grid],[0.5,size(H_variance_norm,2)+0.5],'linewidth',1.5,'color','k')    
end
for y_grid = 0.5:1:size(H_variance_norm,2)+0.5
    line([0.5,size(H_variance_norm,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','k')    
end 

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Permutation_Test_H_occurance',fn_path,1);
    close all;
end

%% Plot the H autocorrelation
figure('position',[680   101   858   877]); hold on; 
H_weight_zscore = zscore(H_weight',0,2);
tau = NaN(1,size(H_weight_zscore,1));
for i = 1:size(H_weight_zscore,1)
    %plot pdf     
    [xc,lags] = xcorr(H_weight_zscore(i,:),20*13,'coeff'); 
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

%%
if savefigs
    handles = get(groot, 'Children');
    saveas(handles,[fn_path 'BinaryState_Vertical'],'svg')
    close all;
end

%% Plot the binary state ID 
figure('position',[462 200 1029 700]); hold on; 
subplot(2,1,1)
imagesc(clusters(:,1:(end-1))',[0 1.25]); 
colormap magma
xlabel('Behavioral State')
set(gca,'YTickLabel',labels,'TickLength',[0,0],'XAxisLocation','top','XTick',[])
setFigureDefaults
set(gca,'position',[8 12 15 4],'box','on')
for x_grid = 0.5:1:size(clusters,1)+0.5
    line([x_grid,x_grid],[0.5,numel(labels)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid =0.5:1:numel(labels)+0.5
    line([0.5,size(clusters,1)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  

subplot(2,1,2)
imagesc(expvar_state,[0 1]); 
colormap magma
xlabel('Behavioral State')
set(gca,'YTick',(1:numel(expvaridx)),'YTickLabel',arrayfun(@(x) num2str(x),(1:numel(expvaridx)),'UniformOutput',0),...
    'TickLength',[0,0])
ylabel('Basis Motifs')
c = colorbar;
set(c,'units','centimeters','position',[24 3 0.5 8.5])
ylabel(c,'Percent Explained Variance')
setFigureDefaults
set(gca,'position',[8 3 15 8.5],'box','on')

%% Just the explained variance
figure; hold on;
imagesc(expvar_state,[0 1]); 
colormap magma
xlabel('Behavioral State')
set(gca,'YTick',(1:numel(expvaridx)),'YTickLabel',arrayfun(@(x) num2str(x),(1:numel(expvaridx)),'UniformOutput',0),...
    'TickLength',[0,0])
ylabel('Basis Motifs')
ylim([0.5,numel(expvaridx)+0.5])
xlim([0.5,size(clusters,1)+0.5])

for x_grid = 1.5:1:size(clusters,1)
    line([x_grid,x_grid],[0.5,numel(expvaridx)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid =1.5:1:size(clusters,1)+0.5
    line([0.5,size(clusters,1)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
c = colorbar;
set(c,'units','centimeters','position',[12.5 3 0.5 8.5])
ylabel(c,'Percent Explained Variance')
setFigureDefaults
set(gca,'position',[3 2 8.5 8.5],'box','on')

%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Motifs vs State',fn_path,1);
    close all;
end

%% save off consecutive example frames (allow up to 5 missed frames) (max 5 instances) for each cluster
frame_indices = arrayfun(@(x) find(indx_clusters==x),unique_states,'UniformOutput',0); 
% get starts and stops of consecutive indices
consecutive_frames = cellfun(@(x) find(diff([false;[1;diff(x)]==1;false])~=0), frame_indices,'UniformOutput',0);
consecutive_frames = cellfun(@(x) reshape(x', 2,[])', consecutive_frames,'UniformOutput',0);
% select up to 5 longest instances
[t, temp] = cellfun(@(x) maxk(x(:,2)-x(:,1),min([10,size(x,1)])),consecutive_frames,'UniformOutput',0);
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





% % % features_reconstructed = NaN(size(features_downsampled));
% % % behav_state_loadings = NaN(1,numel(unique_states));
% % % for i = 1:numel(unique_states)
% % %     temp =  NaN(size(features_downsampled));
% % %     temp(indx_clusters==unique_states(i),:) = repmat(nanmean(features_downsampled(indx_clusters==unique_states(i),:),1),sum(indx_clusters==unique_states(i)),1); 
% % %     features_reconstructed(indx_clusters==unique_states(i),:) = repmat(nanmean(features_downsampled(indx_clusters==unique_states(i),:),1),sum(indx_clusters==unique_states(i)),1);    
% % %     
% % %     temp(isnan(temp))=0;
% % %     behav_state_loadings(i) = 1 - nanvar(features_downsampled(:)-temp(:))./nanvar(features_downsampled(:));
% % % end
% % % 
% % % %zero out any nans since otherwise they won't be factored into the expvar
% % % features_reconstructed(isnan(features_reconstructed))=0;
% % % behav_state_expvar = 1 - nanvar(features_downsampled(:)-features_reconstructed(:))./nanvar(features_downsampled(:));
% % % behav_state_loadings = behav_state_loadings/(sum(behav_state_loadings));
% % % figure('position',[680   458   560   520]); hold on;
% % % bar(behav_state_loadings*100,'facecolor',[0.5 0.5 0.5])
% % % for i = 1:numel(behav_state_loadings)
% % %     text(i,behav_state_loadings(i)*100+0.5,sprintf('%.2g%%',behav_state_loadings(i)*100),'Rotation',90,'FontSize',14,'FontName','Arial','Fontweight','normal')
% % % end
% % % xlim([0.5,numel(unique_states)+0.5])
% % % ylabel({'Relative Percent';'Explained Variance'})
% % % title({'Behavioral State';'Contributions'},'FontName','Arial','FontSize',16,'FontWeight','normal')
% % % xlabel('Behavioral State')
% % % setFigureDefaults;
% % % set(gca,'position',[3 3 6 8.5])
% % % 
% % % %shuffle time-state correlation
% % % rng('default')
% % % behav_state_expvar_shuffled = NaN(1,1000);
% % % for cur_shuf = 1:1000
% % %     features_reconstructed_shuffled = NaN(size(features_downsampled));
% % %     indx_clusters_shuffled = indx_clusters(randperm(numel(indx_clusters)));
% % %     for i = 1:numel(unique_states) %shuffled the labels so applying mean state to incorrect indices
% % %         features_reconstructed_shuffled(indx_clusters==unique_states(i),:) = repmat(nanmean(features_downsampled(indx_clusters_shuffled==unique_states(i),:),1),sum(indx_clusters==unique_states(i)),1);    
% % %     end
% % %     features_reconstructed_shuffled(isnan(features_reconstructed_shuffled))=0;
% % %     behav_state_expvar_shuffled(cur_shuf) = 1 - nanvar(features_downsampled(:)-features_reconstructed_shuffled(:))./nanvar(features_downsampled(:));
% % % end
% % % %right tailed pvalue 
% % % pval = sum([behav_state_expvar_shuffled,behav_state_expvar]>=behav_state_expvar)/numel([behav_state_expvar_shuffled,behav_state_expvar]);
% % % 
% % % figure; hold on;
% % % histogram([behav_state_expvar_shuffled,behav_state_expvar],'numbins',100,'edgecolor','none','facecolor',[0.5 0.5 0.5])
% % % line([behav_state_expvar,behav_state_expvar],[0 max(get(gca,'ylim'))],'linestyle','--','color',[0.75 0 0],'linewidth',2)
% % % title({'Permutation Test with Average';'Clusters from Shuffled Indices'},'Fontweight','normal','FontSize',16,'FontName','Arial');
% % % xlabel('Percent Explained Variance');
% % % ylabel('Shuffle Counts');
% % % setFigureDefaults;
% % % set(gca,'position',[2 2 5 5])
% % % 
% % % %shuffle time-state correlation
% % % rng('default')
% % % behav_state_expvar_shuffled = NaN(1,1000);
% % % for cur_shuf = 1:1000
% % %     features_reconstructed_shuffled = NaN(size(features_downsampled));
% % %     indx_clusters_shuffled = indx_clusters(randperm(numel(indx_clusters)));
% % %     for i = 1:numel(unique_states) %shuffled the labels so applying mean state to incorrect indices
% % %         features_reconstructed_shuffled(indx_clusters_shuffled==unique_states(i),:) = repmat(nanmean(features_downsampled(indx_clusters==unique_states(i),:),1),sum(indx_clusters==unique_states(i)),1);    
% % %     end
% % %     features_reconstructed_shuffled(isnan(features_reconstructed_shuffled))=0;
% % %     behav_state_expvar_shuffled(cur_shuf) = 1 - nanvar(features_downsampled(:)-features_reconstructed_shuffled(:))./nanvar(features_downsampled(:));
% % % end
% % % %right tailed pvalue 
% % % pval = sum([behav_state_expvar_shuffled,behav_state_expvar]>=behav_state_expvar)/numel([behav_state_expvar_shuffled,behav_state_expvar]);
% % % 
% % % figure; hold on;
% % % histogram([behav_state_expvar_shuffled,behav_state_expvar],'numbins',100,'edgecolor','none','facecolor',[0.5 0.5 0.5])
% % % line([behav_state_expvar,behav_state_expvar],[0 max(get(gca,'ylim'))],'linestyle','--','color',[0.75 0 0],'linewidth',2)
% % % title({'Permutation Test with Average';'Clusters Applied to Shuffled Indices'},'Fontweight','normal','FontSize',16,'FontName','Arial');
% % % xlabel('Percent Explained Variance');
% % % ylabel('Shuffle Counts');
% % % setFigureDefaults;
% % % set(gca,'position',[2 2 5 5])
% % % 
% % % %shuffle labels
% % % rng('default')
% % % behav_state_expvar_shuffled = NaN(1,1000);
% % % for cur_shuf = 1:1000
% % %     features_reconstructed_shuffled = NaN(size(features_downsampled));
% % %     unique_states_shuffled = unique_states(randperm(numel(unique_states)));
% % %     for i = 1:numel(unique_states) %shuffled the labels so applying mean state to incorrect indices
% % %         features_reconstructed_shuffled(indx_clusters==unique_states(i),:) = repmat(nanmean(features_downsampled(indx_clusters==unique_states_shuffled(i),:),1),sum(indx_clusters==unique_states(i)),1);    
% % %     end
% % %     features_reconstructed_shuffled(isnan(features_reconstructed_shuffled))=0;
% % %     behav_state_expvar_shuffled(cur_shuf) = 1 - nanvar(features_downsampled(:)-features_reconstructed_shuffled(:))./nanvar(features_downsampled(:));
% % % end
% % % %right tailed pvalue 
% % % pval = sum([behav_state_expvar_shuffled,behav_state_expvar]>=behav_state_expvar)/numel([behav_state_expvar_shuffled,behav_state_expvar]);
% % % figure; hold on;
% % % histogram([behav_state_expvar_shuffled,behav_state_expvar],'numbins',100,'edgecolor','none','facecolor',[0.5 0.5 0.5])
% % % line([behav_state_expvar,behav_state_expvar],[0 max(get(gca,'ylim'))],'linestyle','--','color',[0.75 0 0],'linewidth',2)
% % % title({'Permutation Test Shuffled';'Cluster Labels'},'Fontweight','normal','FontSize',16,'FontName','Arial');
% % % xlabel('Percent Explained Variance');
% % % ylabel('Shuffle Counts');
% % % setFigureDefaults;
% % % set(gca,'position',[2 2 5 5])
% % % 
% % % if savefigs
% % %     handles = get(groot, 'Children');
% % %     saveCurFigs(handles,'-svg','Explained Variance Figures',fn_path,1);
% % %     close all;
% % %     save([fn_path 'BehavioralStateStats'],'clusters','pval','behav_state_loadings','behav_state_expvar');
% % % end



