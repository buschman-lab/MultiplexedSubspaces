%Realtime Motif Classification

fn = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Long_Recordings\431-10-17-2019_1Fitted_block_hemoflag0_1';

%load data
data = load(fn,'bad_pxl','data_test');
bad_pxl = data.bad_pxl;
data = data.data_test;
data_full = NaN(numel(bad_pxl),size(data,2));
data_full(bad_pxl==0,:) = data;

%Get the h weightings for each motif (reconvolved for timing)
expvaridx = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\OrderedByExpVar_Final.mat');
expvaridx = expvaridx.idx;
tempdata = load(fn,'w','H');
w = tempdata.w(:,expvaridx,:); 
h = tempdata.H(expvaridx,:);
h_weight = NaN(size(h))';
for cur_motif = 1:size(h,1)
    h_weight(:,cur_motif) = helper.reconstruct(nanmean(w(:,cur_motif,:),1),h(cur_motif,:)); 
end %motif loop
clear tempdata w
%%
%threshold
% [h_thresh, h_bin, h_onset] = ThresholdMatrix(h_weight',2,'std');
[h_thresh, h_bin, h_onset] = ThresholdMatrix(h,1.5,'std');

%remove onset too close to begining of recording (first 30 seconds and last 30 seconds)
h_onset = cellfun(@(x) x(x>=13*30), h_onset, 'UniformOutput',0);
h_onset = cellfun(@(x) x(x<(size(data,2)-13*30)), h_onset, 'UniformOutput',0);
%%
%get the average image before each motif
indx = (-13:5);
mean_image = reshape(nanmean(data_full,2),[68,68]);
mean_image = (mean_image-nanmin(mean_image(:))/nanmax(mean_image(:))-nanmin(mean_image(:)));
mean_image_pre = NaN(size(data_full,1),numel(indx) ,numel(h_onset));
figure; 
for i = 1:numel(h_onset)     
   for j = 1:numel(indx)      
      mean_image_pre(:,j,i) = nanmean(data_full(:,h_onset{i}+indx(j)),2);             
      temp = reshape(nanmean(data_full(:,h_onset{i}+indx(j)),2),[68 68]);   
%       temp = (temp-nanmin(temp(:)))/(nanmax(temp(:))-nanmin(temp(:)));
      imagesc(temp,[0 0.05]);
      title(sprintf('motif %d : frame %d',i,indx(j)));
      axis square
      pause(0.2)
   end       
end

%% Binary classification. Motif x vs all others
rng('default');
indx = (-13:36);
auc = NaN(numel(h_onset),numel(indx));
feature_pixels = cell(numel(h_onset),numel(indx));
for cur_t = 1:numel(indx) %timepoint loop
   fprintf('\nWorking on timepoint %d',cur_t);
   features = cell(1,numel(h_onset));
   for cur_m = 1:numel(h_onset) %motif loop
       %get features for each motif
       features{cur_m} = data(:,h_onset{cur_m}+indx(cur_t));
   end
   
   %binarize features each motif
   for cur_m = 1:numel(h_onset) %motif loop
      fprintf('\n\t motif %d',cur_m);
      num_features_split = floor(size(features{cur_m},2)/(numel(h_onset)-1)); %evenly sample all subgroups
      num_features = num_features_split*(numel(h_onset)-1);
      %trim first group as necessary so can be evenly split
      group_one = features{cur_m}(:,randperm(size(features{cur_m},2),num_features));
      
      group_two = cell(1,numel(h_onset));
      off_target_indx = find(((1:numel(h_onset))~=cur_m)==1);
      for i = off_target_indx
         group_two{i} = features{i}(:,randperm(size(features{i},2),num_features_split));
      end
      group_two=group_two(~cellfun('isempty',group_two));
      temp = cat(2,group_one,group_two{:})';
      
%       %%PIXELWISE
%       features_parsed = cat(2,temp,cat(1,ones(num_features,1),ones(num_features,1)*2)); 
%       %loop through each pixel       
%       pxl_auc = NaN(1,size(temp,2));
%       for cur_pxl = 1:size(temp,2)          
%            [~, Observed, ~,~] = SVMClassifier_Binary(features_parsed(:,[cur_pxl,end]),[],'holdout',0.2,'verbose',0,...
%                'nshuf',0,'featureselect','none','numberfeatures',50,'optimize',0,'optimize_maxiter',20,...
%                'pca',0,'solver',1,'kernel','rbf','numkfold',5);
%            pxl_auc(cur_pxl) = Observed.AUC;           
%       end
%       temp_px_full = zeros(numel(bad_pxl),1);
%       temp_px_full(bad_pxl==0)=pxl_auc;
%       temp_px_full(temp_px_full==0)=NaN;
%       feature_pixels{cur_m,cur_t} = temp_px_full;

      %%ANOVA
      features_parsed = cat(2,temp,cat(1,ones(num_features,1),ones(num_features,1)*2));                  
      [px, Observed, ~,TrainAUC] = SVMClassifier_Binary(features_parsed,[],'holdout',0.2,'verbose',0,...
       'nshuf',0,'featureselect','anova','numberfeatures',150,'optimize',0,'optimize_maxiter',20,...
       'pca',0,'solver',1,'kernel','rbf','numkfold',5);
      auc(cur_m,cur_t) = Observed.AUC; 
      temp_px = zeros(size(data,1),1);
      temp_px(px)=1;
      temp_px_full = zeros(numel(bad_pxl),1);
      temp_px_full(bad_pxl==0)=temp_px;
      temp_px_full(temp_px_full==0)=NaN;
      feature_pixels{cur_m,cur_t} = temp_px_full;
   end      
   
end

%%
kernel = 3;
auc_smooth = auc';
for i =1:size(auc,1)
   auc_smooth(:,i) = movmean(auc_smooth(:,i),kernel);
end
auc_smooth = auc_smooth';
figure; plot((1:numel(indx)),auc_smooth','linewidth',2);
% figure; plot(auc');
set(gca,'xtick',1:5:numel(indx))
set(gca,'xticklabel',indx(1:5:end)*75,'xticklabelrotation',30)
line([find(indx==0),find(indx==0)],[0 1],'color','r','linewidth',2,'linestyle','--')
ylim([0.4 1])
xlabel('time (ms)')
setFigureDefaults

% [h,p] =ttest(auc,0.5);
% line([find(h==1,1,'first'),find(h==1,1,'last')],[0.95 0.95],'color',[0.5, 0.5, 0.5],'linewidth',2)


%% Get the weigthed image of pixels used. 
timepoints = 10:13;
figure; 
for cur_m = 8
   %get average weighted images       
   avg_weight = reshape(nanmean(cat(2,feature_pixels{cur_m,timepoints}).*cat(2,auc(cur_m,timepoints)),2),[68 68]);
   sum_weight = reshape(nansum(cat(2,feature_pixels{cur_m,timepoints}).*cat(2,auc(cur_m,timepoints)),2),[68 68]);
   cla
   imagesc(sum_weight);
   colorbar
   axis square
   axis off
   colormap magma
   pause
end

%%


%%
for cur_t = 1:11
   for cur_m = 1:3
   temp = feature_pixels{cur_m,cur_t};
   temp(temp==0)=NaN;
   feature_pixels{cur_m,cur_t} = temp;
   end
end

%% make the image for each motif showing it's most used pixels













































%%
%Load data
D = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TestRepitoires\Clustered_Fit_Kval28_Lambda1';
temp = gatherData(318,'Fitted_block_ModFlag0',0,D,{'H'});
%just use from one mouse
mouse = [temp(:).Mouse];
H = {temp(:).H};
H = H(mouse==mouse(1));

motif = 7; 
%look through the recordings and threshold the Hs
trace = {};
for i = 1:numel(H)
    onset = zeros(size(H{i}));
    for cur_H = 1:size(H{i},1)
        temp = H{i}(cur_H,:);
        temp(temp>(nanmean(temp,2)+(1*std(temp,[],2)))) = 1;
        temp(temp~=1) = 0; 
        temp = diff(temp); 
        temp(temp~=1)=0;
        temp = [temp 0];%pad to return same size
        onset(cur_H,:) = H{i}(cur_H,:).*temp; %get the value at that point        
    end %h loop
    bad_idx = sum(onset,1)==0; 
    %Just take the highest values and create a vector of onset times
    [b,indx] = max(onset,[],1);
    indx(bad_idx)=NaN;   
    %Also set the first second to zero since you won't be able to use those
    %timepoints and look retrospectively. 
    indx(1:13)=NaN;
    indx(end-13:end)=NaN;
    trace{i} = indx;
end

% Get the data point before each motif (reconstruct)
tempdata = gatherData(318,'Fitted_block_ModFlag0',0,D,{'bad_pxl','data_test'});
tempdata = tempdata(mouse==mouse(1));
alldata = {};
for i = 1:numel(tempdata)
    temp = zeros(size(tempdata(i).bad_pxl,1),size(tempdata(i).data_test,2));
    temp(~tempdata(i).bad_pxl,:)=tempdata(i).data_test;    
    alldata{i}=temp;
end

%Concatenate alldata and the trace
trace = cat(2,trace{:});
trace(isnan(trace))=0;
alldata = cat(2,alldata{:});


%% Run the classification
Observed = {};
COUNT = 1; 
offsets = [0:1:8];
for nframes=offsets
rng('default')
fprintf('\nWorking on classification round %d, nframes=%d',COUNT,nframes);
%Sweep offsets
begin = 1; %-1 = 2 prior to start (I think)
%Get the vectors prior to each event and the motif it predates
motifs = unique(trace);
motifs(motifs==0)=[];

predictors = {};
label = {};
for i = 1:numel(motifs)
    index = find(trace==motifs(i));
    for cur_ind = 1:numel(index)
        start = index(cur_ind)-begin;
        temp = alldata(:,(start:start+nframes));
%         start = index(cur_ind)+nframes;
%         temp = alldata(:,(start:start+2));
        temp = temp(:);
        predictors{i}(cur_ind,1:numel(temp)) = temp;        
        label{i}(cur_ind)=i;
    end
end
   
%Now create a classification matrix
feature_matrix = cat(1,predictors{:});
label = cat(2,label{:})';

%Remove columns with zeros variance across all instances
bad_dim = nanvar(feature_matrix,[],1)<=eps;
feature_matrix(:,bad_dim)=[];

%Classify
for i = 1:numel(motifs)
    %binarize
    temp = label;
    temp(temp~=i)=0;
    temp(temp==i)=1;   
    temp = temp+1; %must be positive
    [~, Observed{COUNT}(i), ~] = SVMClassifier_Binary_linear([feature_matrix,temp],'holdout',0.25,'featureselect','none','numberfeatures',1000,'pca',0);
end
COUNT = COUNT+1;
end

AUC = [];
for i =1:numel(Observed)
temp = Observed{i};
temp = [temp(:).AUC];
AUC{i}=temp';
end
AUC = cat(2,AUC{:});
AUC = AUC*100;

%%
figure; hold on
x = (offsets(1:size(AUC,2))-1)*75;
y = nanmean(AUC,1);
err = nanstd(AUC,[],1)/sqrt(size(AUC,1));
shadedErrorBar(x,y,err);
line([0 0],[0 100],'linewidth',2,'color','r','linestyle','--');
setFigureDefaults
ylim([50 75]);
ylabel({'Classification';'Accuracy (AUC)'});
xlabel({'Predictor Timepoints';'(ms Respective to Motif Onset)'});
xlim([-100 400])
for i =1:size(AUC,2)
   [h,p] = ttest(AUC(:,i),50);
   if p<=0.05/7
       line([x(i)-75,x(i)],[74 74],'linewidth',2,'color',[0.7 0.7 0.7])
   end
end
set(gca,'units','centimeters','position',[3 3 5 8.5])
set(gcf,'position',[680   486   560   492])

















