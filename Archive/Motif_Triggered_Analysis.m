%add paths
addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\'));
%set filepaths
fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse431_10_17_2019\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_dlc = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse9031_10_17_2019';
fn_widefield = '431-10-17-2019_1Fitted_block_hemoflag0_1';

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

savefigs = 0; 

%% load motif H weightings 
H = load([fn_path, fn_widefield],'H');
H = H.H; 
num_frames = size(H,2);

%Convolve with motif intensity to get correct onset timings
load('AverageDPs_1.mat')
weights = (nanmean(W_clust_smooth,1));
H_weight = NaN(size(H));
for i = 1:size(H,1)
    H_weight(i,:) = helper.reconstruct(weights(:,i,:),H(i,:));
end
H_weight = H_weight(expvaridx,:);
[h_thresh, h_bin, h_onset] = ThresholdMatrix(H_weight,bp.min_std);

save([fn_path, 'H_stats'],'h_onset');
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
    saveCurFigs(gcf,'-svg','Chosen ROIs',fn_path,1);
    close; 

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
% features = cat(2,face_motion_energy{:},tail_to_nose,frontpaws_to_nose,limb_speed);
% labels = {'nose motion energy','whisker motion energy','tail to nose','front paws to nose','front right paw speed',...
% 'front left paw speed','back right paw speed','back left paw speed'};
% labels_abbrev = {'NME','WME','T2N','F2N','FRPS','FLPS','BRPS','BLPS'};
features = cat(2,face_motion_energy{:},tail_to_nose,frontpaws_to_nose,limb_speed);
labels = {'nose motion energy','whisker motion energy','tail to nose','front paws to nose','limb speed'};
labels_abbrev = {'NME','WME','T2N','F2N','LS'};




%Trim to match start and stop of imaging 
features = features(onset:offset,:);

% Downsample to match motif duration
features_downsampled = NaN(num_frames,size(features,2));
for i = 1:size(features_downsampled,2)
    temp = features(:,i);
    features_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_frames),'linear');
end

if bp.zscore %optional zscore
   features_downsampled = zscore(features_downsampled,1);
end

clear raw_data facecam_data W_clust_smooth
%% Plot statistics about the behavioral traces
num_features = size(features_downsampled,2);
col = getColorPalet(num_features);
%distributions
figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
for i = 1:num_features
    subplot(r,c,i); hold on;
    %plot pdf 
    temp = features_downsampled(:,i);
%     if i > num_features-4
%        temp(temp<=5)=NaN;
%     end
    [f,xi] = ksdensity(temp); 
    plot(xi,f,'linewidth',2,'color',col(i,:)); 
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;
    if mod(i,3)==1
        ylabel('Probability')
    end    
    if i > num_features-c
        xlabel('Z-Score')
    end
end


% Plot the autocorrelation of each behavioral factor
figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
for i = 1:num_features
    subplot(r,c,i); hold on;
    %plot pdf 
    [xc,lags] = xcorr(features_downsampled(:,i),120*13,'coeff'); 
    idx = [floor(numel(lags)/2):numel(lags)];
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',col(i,:));   
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;    
    ylim([min(xc) 1]);
    xlim([0 max(lags/13)])
    if mod(i,2)
        ylabel('Rho')
    end
    if i > num_features-c
        xlabel('time (s)')
    else
        set(gca,'Xtick',[])
    end
end

% Correlation between factors
figure; hold on; 
rho = corr(features_downsampled);
imagesc(rho,[0 1]);
colorbar
colormap magma
set(gca,'XTick',(1:num_features),'YTick',(1:num_features),'YTickLabel',labels_abbrev,'XTickLabel',labels_abbrev,'XTickLabelRotation',45)
xlim([0.5,num_features+0.5])
ylim([0.5,num_features+0.5])
axis square
title('Correlation Between Behavioral Features','FontSize',16,'FontName','Arial','FontWeight','normal')
setFigureDefaults; 

%% Manually define states
%define split index 
figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
split_idx = {[-0.8],[-0.5],[-1.75],[-3],[1]};
features_binned = NaN(size(features_downsampled));
for i = 1:num_features
    temp = features_downsampled(:,i);
    temp(temp<split_idx{i})=-100;
    temp(temp>split_idx{i})=100;
    features_binned(:,i) = temp;
    
    subplot(r,c,i); hold on;
    plot(temp)
end

%% Compute the relative explained variance of each motif during each behavioral cluster
[clusters,~,indx_clusters] = unique(features_binned,'rows');
indx_clusters = movmode(indx_clusters,floor(13*2.5));
temp = [];
unique_clust = unique(indx_clusters);
for i = 1:numel(unique_clust)
    temp(i) = sum(indx_clusters==unique_clust(i));  
end
clusters = clusters(unique_clust,:);
clusters(:,end+1) = temp;

%% Calculate Motif-Triggered Snippet Clusters
snippets_clusters = cell(1,size(H_weight,1));
for cur_motif = 1:size(H_weight,1)    
    %remove any motif onsets that are too close to start and end
    temp = h_onset{cur_motif};
    temp = temp(temp-abs(min(bp.trig_dur))>=1 & temp+abs(max(bp.trig_dur))<num_frames);
    temp = arrayfun(@(x) indx_clusters(x+bp.trig_dur), temp, 'UniformOutput',0);       
    snippets_clusters{cur_motif} = cat(2,temp{:});
end

%% For each motif snippet get the most frequent behavioral cluster
temp = cellfun(@(x) (mode(x,1)), snippets_clusters,'UniformOutput',0);
group = cellfun(@(n) ones(1,numel(temp{n}))*n, num2cell(1:size(H_weight,1)),'UniformOutput',0);
temp = [temp{:}];
group = [group{:}];

%unique states
unique_states = unique(temp);
motifs = (1:14);
data = NaN(numel(motifs),numel(unique_states));
for i = 1:numel(motifs)
    for j = 1:numel(unique_states)
        data(i,j) = sum(temp==unique_states(j) & group==motifs(i));
    end    
end

%% Look at the normalized number of triggers related to each state
figure;
temp = data./nansum(data,2);
data_normcol = temp./nansum(temp,1);
imagesc(data_normcol,[0 0.15]);
colormap magma

%% Look at the relative frame-wise explained variance of each motif for each state
data = load([fn_path fn_widefield],'w','data_test','H');
w = data.w(:,expvaridx,:);
H = data.H(expvaridx,:);

load_frame = NaN(size(H,1),num_frames);
for cur_motif = 1:size(H,1)
    fprintf('\nworking on motif %d',cur_motif);
    WH = helper.reconstruct(w(:,cur_motif,:),H(cur_motif,:));  
    for cur_f = 1:num_frames
        load_frame(cur_motif,cur_f) = 1 - nanvar(data.data_test(:,cur_f)-WH(:,cur_f))./nanvar(data.data_test(:,cur_f));        
    end   
end %motif loop
clear data w H WH

%remove motifs that capture less than 1 percent of the variance to a given frame
load_frame(load_frame<0.01)=NaN;

%split by behavioral state
expvar_state = arrayfun(@(x) load_frame(:,indx_clusters==x),unique_states,'UniformOutput',0); 
expvar_state = cellfun(@(x) nansum(x./nansum(x(:)),2), expvar_state, 'UniformOutput',0);
expvar_state = [expvar_state{:}];
        
figure;
imagesc(expvar_state,[0 0.15]);
colormap magma

%%



























%% Calculate Motif-Triggered Snippets



snippets = cell(1,size(H_weight,1));
for cur_motif = 1:size(H_weight,1)    
    %remove any motif onsets that are too close to start and end
    temp = h_onset{cur_motif};
    temp = temp(temp-abs(min(bp.trig_dur))>=1 & temp+abs(max(bp.trig_dur))<num_frames);
    temp = arrayfun(@(x) features_downsampled(x+bp.trig_dur,:), temp, 'UniformOutput',0);       
    snippets{cur_motif} = cat(3,temp{:});
end

alldata = cat(3,snippets{:});
motif_triggered_data = NaN(size(alldata,3)*size(alldata,1),size(alldata,2));
for i = 1:size(alldata,2)
    temp = squeeze(alldata(:,i,:));
    motif_triggered_data(:,i) = temp(:);
end
clear alldata;

%% Plot statistics about the behavioral traces
num_features = size(motif_triggered_data,2);
col = getColorPalet(num_features);
%distributions
figure('position',[680   101   858   877]); hold on; 
[r,c] = numSubplot(num_features,2);
for i = 1:num_features
    subplot(r,c,i); hold on;
    %plot pdf 
    [f,xi] = ksdensity(motif_triggered_data(:,i)); 
    plot(xi,f,'linewidth',2,'color',col(i,:)); 
    title(sprintf('%s',labels{i}),'FontName','Arial','FontSize',16,'FontWeight','normal')
    setFigureDefaults;
    if mod(i,3)==1
        ylabel('Probability')
    end    
    if i > num_features-c
        xlabel('Z-Score')
    end
end




%% Plot average motif triggered, grouped by limb
col = getColorPalet(size(H_weight,1));
col = num2cell(col',1);
for i = 1:size(snippets{1},2)
    figure('position',[315   558   925   420]); hold on; 
    pt = cellfun(@(x,y) Plot_Snippet(squeeze(x(:,i,:)),(bp.trig_dur*75)',y), snippets,col,'UniformOutput',0);
    xlabel('Time (ms)');
    ylabel('Z-Score')
    title(labels{i},'FontName','Arial','FontWeight','normal','Fontsize',16);
    legend([pt{:}],cellfun(@(x) num2str(x),(num2cell(1:size(H_weight,1))),'UniformOutput',0),'location','eastoutside')
    xlim([min(bp.trig_dur*75),max(bp.trig_dur*75)]);   
    setFigureDefaults;    
    set(gca,'position',[2 2 12 8]);
end

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GroupedByLimb',fn_path,1);
    close all;  
end

%% Plot the average motif triggered, grouped by motif
col = getColorPalet(size(snippets{1},2));
col = num2cell(col',1);
for i = 1:size(H_weight,1)
    figure('position',[315   558   925   420]); hold on; 
    pt = arrayfun(@(n) Plot_Snippet(squeeze(snippets{i}(:,n,:)),(bp.trig_dur*75)',col{n}), 1:size(snippets{i},2),'UniformOutput',0);
    xlabel('Time (ms)');
    ylabel('Z-Score') 
    legend([pt{:}],labels,'location','eastoutside')
    title(sprintf('Motif %d',i),'FontName','Arial','FontWeight','normal','Fontsize',16);
    setFigureDefaults;    
    set(gca,'position',[2 2 12 8]);
end
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GroupedByMotif',fn_path,1);
    close all;  
end

%% Plot the average motif triggered, grouped by motif
col = getColorPalet(size(snippets{1},2));
col = num2cell(col',1);
ylist{[-.225,.225],[-.225,.225],[-.3,.3]};
for i =  1:size(H_weight,1)
    figure('units','centimeters','position',[0 0 15 35]); hold on;
    ax =[];
    for j = 1:5 %size(snippets{i},2)
        ax(j)=subplot(size(snippets{i},2),1,j); hold on
        if j~=5
            Plot_Snippet(squeeze(snippets{i}(:,j,:)),(bp.trig_dur*75)',col{j})
        else
            Plot_Snippet(squeeze(snippets{i}(:,5,:)),(bp.trig_dur*75)',col{5})
            Plot_Snippet(squeeze(snippets{i}(:,6,:)),(bp.trig_dur*75)',col{6})
            Plot_Snippet(squeeze(snippets{i}(:,7,:)),(bp.trig_dur*75)',col{7})
            Plot_Snippet(squeeze(snippets{i}(:,8,:)),(bp.trig_dur*75)',col{8})
        end
        line([min(bp.trig_dur*75),max(bp.trig_dur*75)],[0 0],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2);
        xlim([min(bp.trig_dur*75),max(bp.trig_dur*75)]);   
        ylim(ylist{i});
        setFigureDefaults;   
        set(gca,'XColor','none') 
    end
    for j = 1:5
        set(ax(j),'position',[3 2+((j-1)*3.5) 10 3])
    end
end

%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GroupedByMotif_Split',fn_path,1);
    close all; 
end


figure; hold on
plot(time_vec(start_idx+(win_size+1)/2),nanmean(pval_smooth,1),'color',[0.5 0.5 0.5],'linewidth',2,'linestyle','--');
plot(time_vec(start_idx+(win_size+1)/2),fstat_smooth')
plot(time_vec(start_idx+(win_size+1)/2),nanmean(fstat_smooth,1),'k','linewidth',2);

%%











% % 
% % %%  Calculate Motif-Triggered-Motif Snippets
% % snippets_motifs = cell(1,size(H_weight,1));
% % for cur_motif = 1:size(H_weight,1)    
% %     %remove any motif onsets that are too close to start and end
% %     temp = h_onset{cur_motif};
% %     temp = temp(temp-abs(min(bp.trig_dur))>=1 & temp+abs(max(bp.trig_dur))<num_frames);
% %     temp = arrayfun(@(x) H_weight(:,x+bp.trig_dur), temp, 'UniformOutput',0);       
% %     snippets_motifs{cur_motif} = cat(3,temp{:});
% % end
% % 
% % %% Plot motif triggered motif snippets
% % col = getColorPalet(size(snippets_motifs{1},2));
% % col = num2cell(col',1);
% % for i = 1:size(H_weight,1)
% %     figure('position',[315   558   925   420]); hold on; 
% %     pt = arrayfun(@(n) Plot_Snippet(squeeze(snippets_motifs{i}(n,:,:)),(bp.trig_dur*75)',col{n}), 1:size(snippets_motifs{i},1),'UniformOutput',0);
% %     xlabel('Time (ms)');
% %     ylabel('Weight (au)') 
% %     legend([pt{:}],arrayfun(@(x) num2str(x),(1:size(snippets_motifs{1},2)),'UniformOutput',0),'location','eastoutside')
% %     title(sprintf('Motif %d',i),'FontName','Arial','FontWeight','normal','Fontsize',16);
% %     setFigureDefaults;    
% %     set(gca,'position',[2 2 12 8]);
% % end
% % 








