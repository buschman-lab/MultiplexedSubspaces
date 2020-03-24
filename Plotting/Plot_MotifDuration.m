function Plot_MotifDuration(fn)
% Camden MacDowell - timeless
% plot the one second motif duration, all the motif durations, and compares
% data is a 1 x N cell array, where each c

% optionally load data. Hardcoded for 2020 manuscript
if nargin < 1
   addpath('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO')
   mouse_num = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO\MouseNumInfoByRec.mat');
   group = isVPA(mouse_num.mousenum);  
   fn = GrabFiles('block',0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_WindowLengthSweep_FullData'});
   %remove saline
   fn(group==1)=[];    
end

%% Compile Data
data = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DFF_Data\AllData_binned_SmallMask_3x_2minTraining.mat','data');
data = data.data(group==0);
%%
temp = cellfun(@(x) nanmean(x),data,'UniformOutput',0);
thresh = linspace(0.01,1,100);
rawdata_peakwidth_avg = NaN(numel(temp),numel(thresh));
for i = 1:numel(temp)
    for j = 1:numel(thresh)
       [~,~,temp_peakwidth] = findpeaks(temp{i},(1:numel(temp{i}))*75,'WidthReference','halfheight','Threshold',thresh(j)*std(temp{i})); 
       rawdata_peakwidth_avg(i,j) = mean(temp_peakwidth);
    end
end
%% 
figure; hold on; 
fp = fig_params; 
plot(thresh,nanmean(rawdata_peakwidth_avg),'Marker','o','color',fp.c_discovery,'markersize',6,'LineWidth',1.5)
title('Activity Burst Duration');
xlabel('Threshold Level (std of pixel intensity)')
ylabel('Mean Peak Width')
setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);
%%
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\ManuscriptRevisionFigures_currentbio';
save([savedir '\burstduration.mat'],'rawdata_peakwidth_avg','-v7.3');
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','rawdata_burstduration',savedir,1);
close all;


%% plot just the 13 second version using the original motifs
w = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\ClusteredDPs_rawdata.mat','W_all');
w = w.W_all;

weight = NaN(size(w,2),size(w,3));
for j = 1:size(w,2) 
  w_flat = nanmean(squeeze(w(:,j,:))); %normalized intensity (keep zeros so that weights by number of active pixels and intensity of those pixels
  weight(j,:) = w_flat/max(w_flat);      
end       
   
figure; hold on; 
x = (1:1:size(weight,2))*75;
err = nanstd(weight,[],1)./sqrt(sum(~isnan(weight),1));
shadedErrorBar(x,nanmean(weight),err)
title({'Average Motif';'Activity for L=13'},'FontSize',14,'Fontweight','normal');
xlabel('Time (s)')
ylabel({'Motif Intensity';'(Normalized)'})
setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);

handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','originalmotifduration',savedir,1);
close all;
%%
weight = cell(1,8);
peakwidth = cell(1,8);
expvar = cell(1,9);

for cur_file = 1:numel(fn)
   if mod(cur_file,round(0.05*numel(fn)))==0
       fprintf('\t%2.0f%% done...\n', (cur_file/numel(fn)*100));
   end
   temp = load(fn{cur_file},'paramsweep');
   temp = temp.paramsweep;
   for i = 1:numel(temp)-1 %loop through L length
      w = temp(i+1).W;
      amp_weight = NaN(size(w,2),size(w,3));
      pw = NaN(1,size(w,2));
      for j = 1:size(w,2) 
          w_flat = nanmean(squeeze(w(:,j,:))); %normalized intensity (keep zeros so that weights by number of active pixels and intensity of those pixels
          amp_weight(j,:) = w_flat/max(w_flat);
          [~,~,temp_peakwidth] = findpeaks(amp_weight(j,:),(1:numel(amp_weight(j,:)))*75,'WidthReference','halfheight','MinPeakDistance',(numel(amp_weight(j,:))-1.1)*75); %set so that you only find the max peak (by setting the minpeakdistance to the length of the motif.
          if ~isempty(temp_peakwidth) %potentially can have no peak if barely any activity until very start or end. this happens on very few motifs (N~1)
              pw(j) = temp_peakwidth;
          end          
      end       
   
      if cur_file == 1 && i==1
          weight{1,i} = amp_weight; %+1 to skip L=1 
          peakwidth{1,i} = pw; %+1 to skip L=1           
      else
          weight{1,i} = cat(1,weight{1,i}, amp_weight); %+1 to skip L=1 
          peakwidth{1,i} = cat(2,peakwidth{1,i}, pw); %+1 to skip L=1 
      end
   end
   
   for i = 1:numel(temp) %get the expvar (include L = 1 here)
      if cur_file == 1 && i==1
          expvar{1,i} = temp(i).ExpVar_all;           
      else
          expvar{1,i} = cat(2,expvar{1,i},  temp(i).ExpVar_all); %+1 to skip L=1 
      end
   end
   
   %get the average duration of each even in the raw data
end %file loop


%% plot the weight trajectory
figure; hold on;
for i = 1:numel(weight)
    x = (1:1:size(weight{i},2))*75;
    err = nanstd(weight{i},[],1)./sqrt(sum(~isnan(weight{i}),1));
    shadedErrorBar(x,nanmean(weight{i}),err)
    drawnow
end
title('Average Activity as a function of Motif length');
xlabel('Time (s)')
ylabel('Normalized Motif Intensity')
setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);

%% Plot the motif durations
temp = MakeCellsEqual(peakwidth,2,1);
temp = cat(1,temp{:})';
dur = [4,7,10,13,26,39,52,65];

figure; hold on;
fp = fig_params;
shadedErrorBar(dur*75,nanmean(temp),nanstd(temp)./sqrt(sum(~isnan(temp))),'lineProps',{'color',fp.c_discovery,...
    'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);
plot(dur*75,nanmean(temp),'o','color',fp.c_discovery,'markersize',6,'LineWidth',1.5)

%add line showing chosen length
line([13*75 13*75],[0 500],'linewidth',fp.dl_line_width,'linestyle','--',...
    'color',[0.1 0.1 0.1])
ylim([150,250])

%add titles and format
xlabel('Motif Length (ms)');
ylabel({'Average Peak';'Width (ms)'});
set(gca,'xtick',[75, 975, 2925, 4875]);
xlim([75 4875]);

setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);


%% Plot the explained variance
temp = cat(1,expvar{:})'*100;
dur = [1,4,7,10,13,26,39,52,65];

figure; hold on;
fp = fig_params;
shadedErrorBar(dur*75,nanmean(temp),nanstd(temp)./sqrt(sum(~isnan(temp))),'lineProps',{'color',fp.c_discovery,...
    'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);
plot(dur*75,nanmean(temp),'o','color',fp.c_discovery,'markersize',6,'LineWidth',1.5)

%add line showing chosen length
line([13*75 13*75],[0 500],'linewidth',fp.dl_line_width,'linestyle','--',...
    'color',[0.1 0.1 0.1])
ylim([75,90])

%add titles and format
xlabel('Motif Length (ms)');
ylabel({'Percent Explained';'Variance'});
set(gca,'xtick',[75, 975, 2925, 4875]);
xlim([75 4875]);

setFigureDefaults;
set(gca,'position',[3,5,6,6])
set(gcf,'position',[680   400  450   600]);

%%
handles = get(groot, 'Children');
saveCurFigs(handles,'-svg','choiceofL',savedir,1);
close all;



%scratch code
% group = {4*75,7*75,10*75,13*75,26*75,39*75};
% figure; hold on;
% boxplot(temp,group,'symbol','');
% plot(nanmedian(temp),'linewidth',2,'color','r');
% ylim([0,prctile(temp(:,end),99.3)])
    
