%add paths
addpath(genpath('C:\Users\macdo\Documents\GitHub\Widefield_Imaging_Analysis'));
% 
fn_path = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse431_10_17_2019\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_bodycam = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
fn_dlc = 'Mouse431_10_17_2019DLC.csv';
fn_savebase = 'Mouse431_10_17_2019';
fn_widefield = '431-10-17-2019_1Fitted_block_hemoflag0_1';
cd(fn_path);
load('wrkplacedata.mat')

% fn_path = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse432_10_17_2019\';
% fn_facecam = 'Cam_0_20191017-171729.avi';
% fn_bodycam = 'Cam_1_20191017-171729_Mouse432_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
% fn_dlc = 'Mouse432_10_17_2019DLC.csv';
% fn_savebase = 'Mouse432_10_17_2019';
% fn_widefield = '432-10-17-2019_1Fitted_block_hemoflag0_1';
% cd(fn_path);
% load('tempdata.mat');

%load behavioral analysis paramters
bp = behavioral_params; 

expvaridx = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;

data = load([fn_path fn_widefield],'w','H','data_test');
w = data.w(:,expvaridx,:); 
H = data.H(expvaridx,:);
data = data.data_test;
Onset_times = cell(1,size(H,1));
Activity = H;
for j = 1:size(H,1)
   XHat = tensor_convolve(w(:,j,:),H(j,:));
   rho_featame = zeros(1,size(XHat,2));
   for i = 1:size(XHat,2)    
      rho_featame(i) = corr(XHat(:,i),data(:,i));
   end 
   [pks,locs] = findpeaks(rho_featame,'MinPeakHeight',0.4,'MinPeakDistance',1*15);
   Activity(j,:) = zeros(size(H(j,:)));
   Activity(j,locs) = pks;
   Onset_times{j} = locs;
end %motif loop

num_featames = size(H,2);
savefigs = 0; 
reprocess = 0; 
fp = fig_params;
% Do a motif triggered limb activity 
%load dlc traced data
limb = load([fn_path,'data.mat']);
feat = limb.features_downsampled';
%% compare with firing rates
[coef, score, latent, ~, explained, mu] = pca(feat');
score = score';
dur = 45;
data = cell(1,numel(Onset_times));
data_pc = cell(1,numel(Onset_times));
for i = 1:numel(Onset_times)
   temp = Onset_times{i};   
   snippet = NaN(size(feat,1),dur,numel(Onset_times));
   snippet_pc = NaN(size(score,1),dur,numel(Onset_times));
   for j = 1:numel(temp)
       try %avoid the end
          start = floor(temp(j)-15);
          snippet(:,:,j) = feat(:,start:start+dur-1); 
          snippet_pc(:,:,j) = score(:,start:start+dur-1); 
       catch
       end
   end
   data{i} = snippet;
   data_pc{i} = snippet_pc;
end

%plot the average between the two
col = getColorPalet(size(feat,1));
figure('position',[28         558        1808         420]); hold on;
for i = 1:numel(data)
    subplot(1,numel(data),i); hold on; 
    x = 1:size(data{i},2);
    y = nanmean(data{i},3);
    err = nanstd(data{i},[],3)/sqrt(size(data{i},3));
    for j = 1:size(y,1)
        shadedErrorBar(x,y(j,:),err(j,:),'lineprops',{'color',col(j,:)});
    end
%     title(sprintf('%s',labels{i}));
%     ylim([-0.4 1.2])    
end

%PC Space
num_pcs = 2;
col = getColorPalet(num_pcs);
figure; hold on;
for i = 1:numel(data_pc)
    subplot(7,2,i); hold on; 
    x = 1:size(data_pc{i}(1:num_pcs,:,:),2);
    y = nanmean(data_pc{i}(1:num_pcs,:,:),3);
    err = nanstd(data_pc{i}(1:num_pcs,:,:),[],3)/sqrt(size(data_pc{i}(1:num_pcs,:,:),3));
    for j = 1:size(y,1)
        shadedErrorBar(x,y(j,:),err(j,:),'lineprops',{'color',col(j,:)});
    end
    title(sprintf('%d',i));
end







































%% Preprocessing Behavioral Videos 
if reprocess
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
        handles = get(groot, 'Children');
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
end
%% Processing DLC Analyzed Data
%load dlc traced data
[~,~,raw_data] = xlsread([fn_path, fn_dlc]);

%get the speed of the 4 limbs
[limbs, id, ~] = parse_dlc(raw_data,{'featontrightpawcenter','featontleftpawcenter','backrightpawcenter','backleftpawcenter'},[],bp.dlc_epsilon);
limb_speed = cellfun(@(x) [0; mean(abs(diff(limbs(:,strcmp(id,x)),1)),2)],{'featontrightpawcenter','featontleftpawcenter','backrightpawcenter','backleftpawcenter'},'UniformOutput',0); 
limb_speed = [limb_speed{:}];
limb_speed = (mean(limb_speed,2));

%combined all features
% face_motion_energy{2} = -1 * face_motion_energy{2}+max(face_motion_energy{2}(:)); %may need to flip the whisker energy if high whisking actually blurs the camera and make low energy  
features = cat(2,face_motion_energy{:},limb_speed);
labels = {'nose motion energy','whisker motion energy','limb speed'};
labels_abbrev = {'NME','WME','LS'};

%Trim to match start and stop of imaging 
features = features(onset:offset,:);

% for i = 1:size(features,2)
%     features(:,i) = convn(features(:,i),ones(130,1)/130,'same');
% end

% Downsample to match motif duration
features_downsampled = NaN(num_featames,size(features,2));
for i = 1:size(features_downsampled,2)
    temp = features(:,i);
    features_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),num_featames),'linear'); 
end

%store the mapping featom original to downsampled 
x_query_ds = linspace(1,size(features,1),num_featames);

clear raw_data facecam_data W_clust_smooth

%log transform limb_speed;
features_downsampled(:,3) = log(features_downsampled(:,3));
features_downsampled = zscore(features_downsampled,1);


%% GMM
rng('default');
Y = features_downsampled;
options = statset('Display','final','MaxIter',500);
gmm = fitgmdist(Y,2,'options',options,'CovarianceType','full','RegularizationValue',0,'SharedCovariance',false);
[clust_idx,nlogL,P,logpdf,d2] = cluster(gmm,Y);

%smooth the probabilties. e.g. get a running product
kernel = 26;
P_smooth = convn(P,ones(kernel,1)/kernel,'same');

%get the max probabilty for each timepoint
[clust_prob,clust_id] = max(log(P_smooth),[],2);
% 
% % %remove points with low probabilties
bad_id = clust_prob<=log(0.50);
clust_prob(bad_id)=[];
clust_id(bad_id)=[];
Y(bad_id,:)=[];


num_states = numel(unique(clust_id));
col = getColorPalet(num_states);

%histogram of the full distribution and each sub_distribution
figure('position',[680   387   458   604]); hold on; 
for i = 1:size(Y,2)
    figure('position',[680   387   458   604]); hold on;
%    subplot(3,1,i); hold on;
   histogram(Y(:,i),'BinWidth',0.1,'FaceAlpha',0.2,'FaceColor','k','EdgeColor','none')
   for j = 1:numel(unique(clust_id))
       histogram(Y(clust_id==j,i),'BinWidth',0.1,'EdgeColor','none','FaceColor',col(j,:),'FaceAlpha',0.7)
   end  
   xlim([-4 4])
   ylabel('Bins');
   if i ==size(Y,2)       
       xlabel('Z-Score')       
   end
   legend('location','EastOutside')
   setFigureDefaults;
   grid off
   set(gca,'position',[3 3 3 4])
end

%% Plot Behavioral State Statistics
%get instances of each state
temp_instances = SplitIntoConsecutiveIndices(clust_id);

%remove any instances shorter than threshold
temp = cellfun(@(x) numel(x)>=5, temp_instances,'UniformOutput',0); %set 1 one to remove empty 
temp = cell2mat(temp);
instances = cell(num_states,1);
for i = 1:num_states
    instances{i} = temp_instances(i,temp(i,:)==1);   
end

temp_median_dur = NaN(3,num_states); temp_num_instances = NaN(1,num_states); temp_timespent = NaN(1,num_states);
for i = 1:num_states   
    temp_timespent(i) = (sum(clust_id == i)/numel(clust_id))*100;       
    temp_num_instances(i) = numel(instances{i});
    temp = cellfun(@ numel, instances{i},'UniformOutput',0);       
    temp_median_dur(1,i) = nanmedian([temp{:}]);
    temp_median_dur(2:3,i) = bootci(1000,@nanmedian,[temp{:}]);    
end

figure; hold on; 
bar(temp_median_dur(1,:),'FaceColor','k','FaceAlpha',0.1,'EdgeColor','k','linewidth',1.5);
errorbar(1:num_states,temp_median_dur(1,:),temp_median_dur(1,:)-temp_median_dur(2,:),temp_median_dur(3,:)-temp_median_dur(1,:),'linewidth',2,'LineStyle','none','color','k')
ylabel({'Number of';'Timepoints'});
title({'Duration of Instances'},'FontName','Arial','FontWeight','normal','FontSize',16);
xlabel('Behavioral State')
setFigureDefaults
set(gca,'position',[4 2 5 6])

figure; hold on; 
bar(temp_timespent,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','k','linewidth',1.5);
title({'Percent';'Total Time'},'FontName','Arial','FontWeight','normal','FontSize',16);
xlabel('Behavioral State')
ylabel('Percent of Timepoints')
setFigureDefaults
set(gca,'position',[4 2 5 6])

figure; hold on; 
bar(temp_num_instances,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','k','linewidth',1.5);
title({'Number of';'Instances'},'FontName','Arial','FontWeight','normal','FontSize',16);
xlabel('Behavioral State')
ylabel('Instances')
setFigureDefaults
set(gca,'position',[4 2 5 6])

stats.state_median_dur = round(temp_median_dur(1,:)/13,2);
stats.state_median_dur_ci = round(temp_median_dur(2:3,:)/13,2);
stats.state_instances = temp_num_instances;
stats.state_timespent = temp_timespent;

%% Get motif weight by state
weight = cell(size(H_weight,2),size(instances,1));
for cur_motif = 1:size(H_weight,2)
   for cur_state = 1:size(instances,1)
      temp = cellfun(@(x) nanmean(H_weight(x,cur_motif)), instances{cur_state},'UniformOutput',0);
      weight{cur_motif,cur_state} = cell2mat(temp);
   end          
end
%flatten and pad 
weight = reshape(weight,size(weight,1)*size(weight,2),1);
weight = MakeCellsEqual(weight,2,1);
weight = reshape(weight,size(H_weight,2),size(instances,1));

%% Compare Weights Across Motifs
%looking at the previous figures, make sure inactive/active in correct row;
% temp = weight;
% weight(:,1) = temp(:,2);%active
% weight(:,2) = temp(:,1);%inactive

fh = Plot_CompareBehaviorStates_2States(weight,[0.8 1.2]);

% save([fn_path,'processed_weights.mat'],'weight');

%% load the processed weights featom each animal
weight =  weight_431;
for i = 1:size(weight,1)
   for j = 1:size(weight,2)
       weight{i,j} = [weight_431{i,j},weight_432{i,j}];
   end       
end
avg_weight = [nanmean(([weight{:,1}])),nanmean(([weight{:,2}]))];
weight_norm = weight;
for i = 1:size(weight,1)
   for j = 1:size(weight,2)
       weight_norm{i,j} = (weight{i,j})/avg_weight(j);
   end       
end

    
fh = Plot_CompareBehaviorStates_2States(weight_norm,[0.9 1.1]);
fh = Plot_CompareBehaviorStates_2States(weight,[0.9 1.1]);

% x = cell2mat(weight(:,1))';
% y = cell2mat(weight(:,2))';
% stats = Plot_CompareMotifStatistic(x,y,'yvals',[0, 7e-4],'xvals',[0, 7e-4]);

    
%% Classify
rng('default');
num_xval = 100; 
auc = nan(1,num_xval);
auc_train = zeros(1,num_xval);

for i = 1:num_xval
    x = cell2mat(weight(:,1));
    y = cell2mat(weight(:,2));

    %remove row with nans
    x(:,~any(~isnan(x)))=[];
    y(:,~any(~isnan(y)))=[];

    len_min = min([size(x,2),size(y,2)]);

    %grab equal random samples featom both; 
    x = x(:,randperm(size(x,2),len_min));
    y = y(:,randperm(size(y,2),len_min));

    %Assign labels
    temp_resp = NaN(len_min*2,1);
    temp_resp(1:len_min)= 1;   
    temp_resp(isnan(temp_resp(:,1)))=2;

    [~, Observed, ~, trainingAUC] = SVMClassifier_Binary([cat(1,x',y'),temp_resp],[],'holdout',0.2,'nshuf',0,'featureselect','none',...
        'optimize',1,'pca',0,'solver',1,'kernel','rbf','numkfold',5,'optimize_maxiter',75);
    auc(i) = Observed.AUC;
    auc_train(i) = nanmean([trainingAUC(:).AUC]);
end

%% Plot the distribution of classifier performance
figure('position',[680   382   974   596]); hold on; 
histogram(auc,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.4 0.4 0.4],'linewidth',1,'BinWidth',.05); 
ylabel({'Number of Cross-';'Validation Runs'})
xlabel({'Classifier AUC'});
line([0.5 0.5],[0 25],'color','r','linestyle','--','linewidth',1.5)
% line([nanmean(auc),nanmean(auc)],[0 25],'color','r','linestyle','--','linewidth',1.5)
text(nanmean(auc)+0.3,300,['\mu=', num2str(round(nanmean(auc),2))],'FontSize',16,'FontWeight','normal','FontName','Arial');
xlim([0 1])
% ylim([0 25])
title('Motif featequency','FontName','Arial','FontWeight','normal','FontSize',16);
setFigureDefaults
set(gca,'position',[4 3 3.5 8.5])
legend off

stats.auc_over_chance= sum(auc>0.5)/numel(auc)*100;
stats.auc_mean = nanmean(auc)*100;
stats.auc_sem = sem(auc,2)*100;
stats.auc_ci_mean = bootci(1000,@nanmean,auc)*100;
stats.auc_median = nanmedian(auc)*100;
stats.auc_ci_median = bootci(1000,@nanmedian,auc)*100;
stats.auc_all = auc; 

if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GMM_Figsnew','C:\Users\macdo\OneDrive\Buschman Lab\AnalysisCode_Repository\Mesoscale Network Dynamics 2019 Analyses\Behavioral_Classification\',1);
    close all;
    save(['C:\Users\macdo\OneDrive\Buschman Lab\AnalysisCode_Repository\Mesoscale Network Dynamics 2019 Analyses\Behavioral_Classification\' 'stats.mat'],'stats');
end

%% Plot the autocorrelation of behavioral variables
col = getColorPalet(3);
figure('position',[1105 524 433 454]); hold on;
for i = 1:size(Y,2)
    %plot pdf 
    [xc,lags] = xcorr(Y(:,i)-nanmean(Y(:,i)),120*13,'coeff'); 
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

%% Plot autocorrelation of motif temporal weightings
figure('position',[209, 101, 1329, 877]); hold on; 
H_temp = (H_weight-nanmean(H_weight,1))';
tau = NaN(1,size(H_temp,1));
for i = 1:size(H_temp,1)
    [xc,lags] = xcorr(H_temp(i,:),20*13,'coeff');         
    idx = [ceil(numel(lags)/2)+1:numel(lags)];
    plot(lags(idx)/13,xc(idx),'linewidth',2,'color',[0.7 0.7 0.7]);        
    xlim([0 max(lags/13)])
    ylabel('Rho')
    xlabel('time (s)')    
    tau(i)=find(xc(idx)>=0.5*xc(idx(1)),1,'last')/13;    
end
title({'Motif Weighting';'Autocorrelation'},'FontName','Arial','FontSize',12,'FontWeight','normal')    
line([0 max(lags/13)],[0 0],'linestyle','--','linewidth',2,'color',[0.25 0.25 0.25]); 
text(3,0.5,sprintf('\\tau = %.2g +/- %.2gs',nanmean(tau(i)),sem(tau,2)),'FontSize',16,'FontName','Arial')
setFigureDefaults; 
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 6 6]) 
ylim([-0.25 1]);    


%%
if savefigs
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','GMM_FigsNew',fn_path,1);
    close all;
    save([fn_path 'stats.mat'],'stats');
end


%% distribution of motifs per state
% for i =1:num_states
%     temp = (cell2mat(weight(:,i)));      
%     figure('position',[680, 558, 1031, 420]); hold on;     
%     CompareViolins(temp,fp);
%     title(sprintf('Motif Weights During State %d',i),'FontName','Arial','Fontsize',16,'Fontweight','normal')
%     ylim([-20 0])
%     xlabel('Basis Motifs');
%     ylabel('Average Weight (log)');    
%     setFigureDefaults
%     set(gca,'position',[3 3 14 6])    
% end
% 
% p = NaN(1,size(weight,1));
% for i =1:size(weight,1)
%     temp = cat(1,weight{i,:});    
%     [p(i),tbl] = anova1(temp',[],'off');
%     
%     figure(); hold on;     
%     CompareViolins(temp,fp);
%     title(sprintf('Motif %d Weights Across States',i),'FontName','Arial','Fontsize',16,'Fontweight','normal')
%     xlabel('States');
%     ylabel('Average Weight (log)');    
%     setFigureDefaults
% %     set(gca,'position',[3 3 14 6])    
% end
% open p
% 
% for i =1:num_states
%     temp = (cell2mat(weight(:,i)));      
%     figure('position',[680, 558, 1031, 420]); hold on;     
%     CompareViolins(temp,fp);
%     title(sprintf('Motif Weights During State %d',i),'FontName','Arial','Fontsize',16,'Fontweight','normal')
%     ylim([-20 0])
%     xlabel('Basis Motifs');
%     ylabel('Average Weight (log)');    
%     setFigureDefaults
%     set(gca,'position',[3 3 14 6])    
% end
%
% %% fit a GMM to the data
% % [~,Y] = pca(features_downsampled); 
% rng('default'); %for reproducibility
% Y = features_downsampled; 
% max_comp = 4; 
% 
% %set options
% options = statset('Display','final','MaxIter',500);
% 
% %model selection 
% AIC = zeros(1,max_comp);
% BIC = zeros(1,max_comp);
% PEV = zeros(size(Y,2),max_comp);
% GMModel = cell(1,max_comp);
% for i = 1:max_comp
%     GMModel{i} = fitgmdist(Y,i,'options',options,'CovarianceType','full','RegularizationValue',0,'SharedCovariance',false);
%     [clust_idx] = cluster(GMModel{i},Y);
%     AIC(i) = GMModel{i}.AIC; %get the AIC
%     BIC(i) = GMModel{i}.BIC; %get the AIC
%     PEV(:,i) = TrialPEV(Y,clust_idx);    
%     fprintf('AIC: %.4d  avgPEV: %d numcomponents: %d ',AIC(i),nanmean(PEV(i)),i);      
% end
% %normalize the AIC;
% AIC_norm = (AIC-min(AIC))/(max(AIC)-min(AIC));
% BIC_norm = (BIC-min(AIC))/(max(BIC)-min(BIC));
% 
% %save off model fitting statistic
% save([fn_path 'gmm_modelselection.mat'],'AIC','PEV','GMModel');
% 
% %% plot the AIC and PEV; 
% figure; hold on
% yyaxis right
% plot(1:max_comp,nanmean(PEV),'marker','o','linewidth',2,'color',[0.3 0.7 1])
% ylabel({'Percent Explained';'Variance'});
% set(gca,'ycolor','k')
% yyaxis left
% plot(1:max_comp,AIC_norm,'marker','o','linewidth',2,'color',[0.1050 0.319, 0.06]);
% plot(1:max_comp,BIC_norm,'marker','o','linewidth',2,'color',[0.85 0.3723 0.0008],'linestyle','-');
% set(gca,'ycolor','k')
% ylabel({'AIC & BIC';'(normalized)'})
% title('GM Model Selection','FontName','Arial','FontSize',16,'FontWeight','normal')
% setFigureDefaults
% line([8 8],[0 1],'linestyle','--','color','k','linewidth',2)
% set(gca,'position',[4 2 6 6])
% xlim([1,max_comp]);
% ylim([0 1])

%%%% Plot the clustering of the states
% [~,Y_plot] = pca(features_downsampled); 
% Y_plot = run_umap([clust_idx,features_downsampled],'n_neighbors',30,'min_dist',0.2,'label_column',1);
% figure; hold on; 
% for i = 1:num_states
%     idx = clust_idx==i;
%     plot(Y_plot(idx,1),Y_plot(idx,2),'.','color',col(i,:),'markersize',0.01)
% end
% title({'Behavior Measures identify';'Spontaneous Behavioral States'},'FontSize',16,'FontWeight','normal','FontName','Arial')
% set(gca,'TickLength',[0 0])
% setFigureDefaults 
% set(gca,'position',[2 2 8 6])
% grid on
% box on

