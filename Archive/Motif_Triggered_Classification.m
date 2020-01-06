%add paths
addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
addpath(genpath('C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\'));
%set filepaths
fn_path = 'C:\Users\macdo\OneDrive\Buschman Lab\Scratch Data\Mouse494_10_17_2019\';
fn_facecam = 'Cam_0_20191017-184642.avi';
fn_dlc = 'Cam_1_20191017-184642_Mouse494_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse494_10_17_2019';
fn_widefield = '494-10-17-2019_1Fitted_block_hemoflag0_1';

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

savefigs = 1; 
cd(fn_path)
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
limb_speed = cellfun(@(x) [0; sum(abs(diff(limbs(:,strcmp(id,x)),1)),2)],{'frontrightpawcenter','frontleftpawcenter','backrightpawcenter','backleftpawcenter'},'UniformOutput',0);

%get the distance between nose and front paws
[frontpaws_to_nose, ~, ~] = parse_dlc(raw_data,{'frontrightpawcenter','frontleftpawcenter'},'nosetip',bp.dlc_epsilon);
frontpaws_to_nose = sum(abs(frontpaws_to_nose),2);

%get the distance from nose to tail
[tail_to_nose, ~, ~] = parse_dlc(raw_data,{'tailroot'},'nosetip',bp.dlc_epsilon);
tail_to_nose = sum(abs(tail_to_nose),2);

%combined all features
features = cat(2,face_motion_energy{:},tail_to_nose,frontpaws_to_nose,limb_speed{:});
labels = {'nose motion energy','whisker motion energy','tail to nose','front paws to nose','front right paw speed',...
'front left paw speed','back right paw speed','back left paw speed'};

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


%% Calculate Motif-Triggered Snippets
snippets = cell(1,size(H_weight,1));
for cur_motif = 1:size(H_weight,1)    
    %remove any motif onsets that are too close to start and end
    temp = h_onset{cur_motif};
    temp = temp(temp-(60*13)>=1 & temp+(60*13)<num_frames);
    temp = arrayfun(@(x) features_downsampled(x+bp.trig_dur,:), temp, 'UniformOutput',0);       
    snippets{cur_motif} = cat(3,temp{:});
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
    xlim([min(bp.trig_dur(bp.plotting_idx))*75,max(bp.trig_dur(bp.plotting_idx))*75])   
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
    xlim([min(bp.trig_dur(bp.plotting_idx))*75,max(bp.trig_dur(bp.plotting_idx))*75])
    setFigureDefaults;    
    set(gca,'position',[2 2 12 8]);
end
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GroupedByMotif',fn_path,1);
    close all;  
end
%%
col = getColorPalet(size(snippets{1},2));
col = num2cell(col',1);
for i =  1:size(H_weight,1)
    figure('units','centimeters','position',[0 0 15 35]); hold on;
    ax =[];
    ylist = [];
    for j = 1:5 %size(snippets{i},2)
        ax(j)=subplot(size(snippets{i},2),1,j); hold on
        if j~=5
            Plot_Snippet(squeeze(snippets{i}(:,j,:)),(bp.trig_dur*75)',col{j})            
        else
            title(sprintf('Motif %d',i),'FontName','Arial','FontWeight','normal','FontSize',16);
            Plot_Snippet(squeeze(snippets{i}(:,5,:)),(bp.trig_dur*75)',col{5})
            Plot_Snippet(squeeze(snippets{i}(:,6,:)),(bp.trig_dur*75)',col{6})
            Plot_Snippet(squeeze(snippets{i}(:,7,:)),(bp.trig_dur*75)',col{7})
            Plot_Snippet(squeeze(snippets{i}(:,8,:)),(bp.trig_dur*75)',col{8})
        end
        line([min(bp.trig_dur*75),max(bp.trig_dur*75)],[0 0],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2);
        xlim([min(bp.trig_dur(bp.plotting_idx))*75,max(bp.trig_dur(bp.plotting_idx))*75])   
        ylist(j,:) = abs(get(gca,'ylim'));
        setFigureDefaults;   
        set(gca,'XColor','none') 
        
    end
    for j = 1:5
        set(ax(j),'position',[3 2+((j-1)*3.5) 7 3])
        set(ax(j),'ylim',[-1*max(ylist(:)),max(ylist(:))])
    end
end

%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GroupedByMotif_Split',fn_path,1);
    close all; 
end


%% Generate a noise snippets using random sampled intervals
rng('default')
snippets_noise = cell(1,size(H_weight,1));
for cur_motif = 1:size(H_weight,1)    
    %remove any motif onsets that are too close to start and end
    temp = h_onset{cur_motif};
%     temp = temp(temp-abs(min(bp.trig_dur))>=1 & temp+abs(max(bp.trig_dur))<num_frames); %use the same overall indices as the regular
%     temp = arrayfun(@(x) features_downsampled(bp.noise_dur(1),:), temp, 'UniformOutput',0); 
    temp = temp(temp-(60*13)>=1 & temp+(60*13)<num_frames);
    temp = arrayfun(@(x) features_downsampled(x+bp.noise_dur,:), temp, 'UniformOutput',0);       
    snippets_noise{cur_motif} = cat(3,temp{:});
end
combined_data_noise = cat(3,snippets_noise{:});

%% Classify between motifs
group = arrayfun(@(n) ones(1,size(snippets{n},3))*n, (1:numel(snippets)), 'UniformOutput',0);
combined_data = cat(3,snippets{:});
combined_data_noise =  cat(3,snippets_noise{:});

group = [group{:}];
%trim to windowed time period
combined_data = combined_data(bp.classification_idx,:,:);
combined_data_noise = combined_data_noise(bp.classification_idx,:,:);
time_vec = bp.trig_dur(bp.classification_idx)*75;

win_size = 3; %actually window size is win_size+1
start_idx = (1:2:size(combined_data,1)-win_size);

%get a list of all the motif comparisons
[p,q] = meshgrid(1:size(H_weight,1), 1:size(H_weight,1));
p = tril(p-1);
q = tril(q,-1);
pairs = [p(p~=0) q(q~=0)];

rng('default') %for reproducibility
auc = NaN(size(pairs,1),numel(start_idx));
auc_noise = NaN(size(pairs,1),numel(start_idx));
for time_bin = 1:numel(start_idx)
    fprintf('\n Working on time bin %d of %d',time_bin, numel(start_idx));
    temp = squeeze(nanmean(combined_data(start_idx(time_bin):start_idx(time_bin)+win_size,:,:),1));
    temp_noise = squeeze(nanmean(combined_data_noise(start_idx(time_bin):start_idx(time_bin)+win_size,:,:),1));
    for cur_pair = 1:size(pairs,1)  
%         disp(cur_pair)
        x = temp(:,(group==pairs(cur_pair,1)));
        y = temp(:,(group==pairs(cur_pair,2)));
        
        %get equal samples from both
        num = min(size(x,2),size(y,2));       
        x = x(:,randperm(size(x,2),num));
        y = y(:,randperm(size(y,2),num));        
        label = cat(1,ones(num,1),2*ones(num,1))';
        
        [~, observed, shuffled] = SVMClassifier_Binary(cat(1,cat(2,x,y),label)','featureselect','none',...
        'nshuf',0,'kernel','linear','optimize',0,'verbose',0);
        
        auc(cur_pair,time_bin) = observed.AUC;        
              
        %get the noise distribution
        x = temp_noise(:,(group==pairs(cur_pair,1)));
        y = temp_noise(:,(group==pairs(cur_pair,2)));
        
        %get equal samples from both
        num = min(size(x,2),size(y,2));       
        x = x(:,randperm(size(x,2),num));
        y = y(:,randperm(size(y,2),num));        
        label = cat(1,ones(num,1),2*ones(num,1))';        
        
        [~, observed, ~] = SVMClassifier_Binary(cat(1,cat(2,x,y),label)','featureselect','none',...
        'nshuf',0,'kernel','linear','optimize',0,'verbose',0);
        auc_noise(cur_pair,time_bin) = observed.AUC;
    end
end

%% Plot the classification accuracy for each motif vs all the others
%smooth the auc trace
auc_average = NaN(size(H_weight,1),size(H_weight,1)-1);
auc_average_noise =NaN(size(H_weight,1),size(H_weight,1)-1);
comparison_id = NaN(size(H_weight,1),size(H_weight,1)-1);
for i = 1:size(H_weight,1)
    figure; hold on; 
    idx = (pairs(:,1)==i | pairs(:,2)==i);
    y = auc(idx,:); 
    
    y_noise = auc_noise(idx,:);
    plot(time_vec(start_idx),nanmean(y_noise),'LineWidth',2,'color',[0.1 0.1 0.1],'linestyle','--');
    
    line([time_vec(start_idx(1)), time_vec(start_idx(length(start_idx)))],[0.5 0.5],'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2)
    plot(time_vec(start_idx),y');
    plot(time_vec(start_idx),nanmean(y),'LineWidth',2,'color','k');
    xlabel('Time (ms)');
    ylabel('AUC')
    setFigureDefaults;    
    title(sprintf('Motif %d',i),'FontName','Arial','FontWeight','normal','Fontsize',16);   
    auc_average(i,:) = nanmean(y,2);
    auc_average_noise(i,:) = nanmean(y_noise,2);
    
    %get the associated motif in the comparison    
    temp = pairs(idx,:);
    comparison_id(i,:) = temp(temp~=i);  
end
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(gcf,'-svg','ClassificationAccuracy_ByMotif',fn_path,1);
    close all; 
end

%% Plot the average AUC for each motif
col = getColorPalet(size(H_weight,1));
figure('position',[241   426   999   552]); hold on; 

% line([0, size(H_weight,1)+1],[0.5 0.5],'linestyle','--','linewidth',2,'color',[0.5 0.5 0.5])
plot(nanmean(auc_average_noise,2),'color',[0.25 0.25 0.25],'linestyle','--','linewidth',2);

h = NaN(size(H_weight,1),1);
pval = NaN(size(H_weight,1),1);
for i = 1:size(H_weight,1)
   rng('default')
   x = (rand(numel(auc_average(i,:)),1)-0.5)/5; 
   scatter(x+i, auc_average(i,:), 25, col(comparison_id(i,:),:), 'filled')
   [pval(i),h(i)] = signrank(auc_average(i,:),auc_average_noise(i,:),'tail','right');
   AddSig(h(i),pval(i),[i,i,max(auc_average(i,:)),max(auc_average(i,:))],14,0.035,1,60)
end
boxplot(auc_average','color','k','Symbol','k+');
set(findobj(gca,'type','line'),'linew',1,'color',[0 0 0 0.5])
t = title('Motifs Reflect Cortical Activity During Specific Spontaneous Behaviors',...
    'Fontsize',18,'Fontweight','normal','FontName','Arial','position',[7.5 0.745]);
xlabel('Motif');
ylabel('AUC');
ylim([0.40 0.65])

setFigureDefaults;    
set(gca,'position',[3 3 20, 8.5])
%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(gcf,'-svg','AUC_Boxplot',fn_path,1);
    close all; 
end

%% Anova between motifs
group = arrayfun(@(n) ones(1,size(snippets{n},3))*n, (1:numel(snippets)), 'UniformOutput',0);
combined_data = cat(3,snippets{:});
for i = 1:size(combined_data,2)
    for j = 1:size(combined_data,3)
        combined_data(:,i,j) = convn(combined_data(:,i,j),ones(6,1),'same');
    end
end

group = [group{:}];
%trim to windowed time period
combined_data = combined_data(bp.classification_idx,:,:);
time_vec = bp.trig_dur(bp.classification_idx)*75;

win_size = 0; %actually window size is win_size+1
start_idx = (1:1:size(combined_data,1)-win_size);

pval = NaN(size(snippets{1},2),numel(start_idx));
for time_bin = 1:numel(start_idx)
    fprintf('\n Working on time bin %d of %d',time_bin, numel(start_idx));
    temp = squeeze(nanmean(combined_data(start_idx(time_bin):start_idx(time_bin)+win_size,:,:),1))';
    for cur_limb = 1:size(snippets{1},2)
        test = temp(:,cur_limb);
        p = anova1(test,group,'off');
        pval(cur_limb,time_bin) = -log10(p);
    end
end

%%

   





























%%

% % %% If you want to use spock 
% % %use spock
% % spock_save_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Temp_Data\';
% % username = input(' Spock Username: ', 's');
% % password = passcode();
% % s_conn = ssh2_config('spock.princeton.edu',username,password);
% % 
% % group = arrayfun(@(n) ones(1,size(snippets{n},3))*n, (1:numel(snippets)), 'UniformOutput',0);
% % combined_data = cat(3,snippets{:});
% % for i = 1:size(combined_data,2)
% %     for j = 1:size(combined_data,3)
% %         combined_data(:,i,j) = convn(combined_data(:,i,j),ones(13,1)/5,'same');
% %     end
% % end
% % group = [group{:}];
% % %trim to windowed time period
% % combined_data = combined_data(30:79-20,:,:);
% % time_vec = bp.trig_dur(30:79-19);
% % 
% % win_size = 4;
% % start_idx = (1:4:size(combined_data,1)-win_size);
% % 
% % %get a list of all the motif comparisons
% % [p,q] = meshgrid(1:size(H_weight,1), 1:size(H_weight,1));
% % p = tril(p-1);
% % q = tril(q,-1);
% % pairs = [p(p~=0) q(q~=0)];
% % 
% % rng('default') %for reproducibility
% % auc = NaN(size(pairs,1),numel(start_idx));
% % for time_bin = 1:numel(start_idx)
% %     fprintf('\n Working on time bin %d of %d\n',time_bin, numel(start_idx));
% %     temp = squeeze(nanmean(combined_data(start_idx(time_bin):start_idx(time_bin)+win_size,:,:),1));
% %     for cur_pair = 1:14 %size(pairs,1)        
% %         x = temp(:,(group==pairs(cur_pair,1)));
% %         y = temp(:,(group==pairs(cur_pair,2)));
% %         
% %         %get equal samples from both
% %         num = min(size(x,2),size(y,2));       
% %         x = x(:,randperm(size(x,2),num));
% %         y = y(:,randperm(size(y,2),num));        
% %         label = cat(1,ones(num,1),2*ones(num,1))';
% %         
% %         savedata = cat(1,cat(2,x,y),label)';
% %         filename = sprintf('behav_class_pair_%d_timebin_%d.mat',cur_pair,time_bin);
% %         save([spock_save_dir filename],'savedata');
% % 
% %         script_name = WriteBashScript(sprintf('%d_%d',time_bin,cur_pair),'Spock_ClassifyBehavioralWaveForms',{filename},{"'%s'"});
% %         response = ssh2_command(s_conn,...
% %             ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
% %             sprintf('sbatch %s',script_name)]); 
% %         
% %     end
% % end
% % 
% % % Close the connection
% % ssh2_close(s_conn);
% % 
% % %clear personal information
% % clear s_conn password username

% % %% Make Classification Figures
% % cd(spock_save_dir)
% % auc = NaN(size(pairs,1),numel(start_idx));
% % for time_bin = 1:numel(start_idx)
% %     for cur_pair = 1:size(pairs,1)  
% %         try
% %         load([spock_save_dir sprintf('Processed_behav_class_pair_%d_timebin_%d.mat',cur_pair,time_bin)]);
% %         auc(cur_pair,time_bin) = Observed.AUC;
% %         catch
% %             auc(cur_pair,time_bin) = NaN;
% %         end
% %     end
% % end












