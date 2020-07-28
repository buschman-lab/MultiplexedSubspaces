

addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Ephys_Processing');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\FunctionsForLFP');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Plotting');
addpath(genpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Toolboxes\TDTSDK'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));

%load previous motifs
motifs = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\AverageDPs_1.mat','W_clust_smooth');
motifs = motifs.W_clust_smooth;

%load new data; 
[data_norm,nanpxs,data_train, data_test] = ProcessAndSplitData('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_Preprocessed\501_11_21_2019_1dff_combined.mat',[]);
data = conditionDffMat(data_norm',nanpxs,[],[68 38,size(data_norm,2)]);
data = SpatialGaussian(data,[3 3]);

%reformat to the same size at motifs
temp = NaN(68,68,size(data,3));
temp(:,1:38,:)=data;
data = temp; clear temp;
data_flat = reshape(data,68*68,size(data,3));
data_flat(isnan(data_flat))=0;

%get only the pixels that are non masked in each
bad_pxl = (nanvar(data_flat,[],2)<=eps) + (nanvar(reshape(motifs,68*68,14*26),[],2)<=eps) >= 1;
data_flat(bad_pxl,:) = [];
motifs(bad_pxl,:,:) = [];

%fit motifs to the new data
chunk = floor(size(data_flat,2)/10);
H = [];
allstats = {};
for i = 1:10
    if i == 10
    [~,H_temp,stats_fit] = fpCNMF(data_flat(:,((i-1)*chunk)+1:end),'non_penalized_iter',...
        0,'penalized_iter',100,...
        'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
        'ortho_H',1,'W',motifs);        
    else
    [~,H_temp,stats_fit] = fpCNMF(data_flat(:,((i-1)*chunk)+1:i*chunk),'non_penalized_iter',...
        0,'penalized_iter',100,...
        'speed','fast','verbose',1,'lambda',0,'w_update_iter',0,...
        'ortho_H',1,'W',motifs);
    end
    H = [H,H_temp];
    allstats{i} = stats_fit;
end
H_orig = H;


%% IF RELOADING DATA
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Ephys_Processing');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\FunctionsForLFP');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Plotting');
addpath(genpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Toolboxes\TDTSDK'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));

motifs = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_Lambda4e-4_Kval28\AverageDPs_1.mat','W_clust_smooth');
motifs = motifs.W_clust_smooth;

load('PreviousMotifFit.mat','H_orig','data_norm','nanpxs')

data = conditionDffMat(data_norm',nanpxs,[],[68 38,size(data_norm,2)]);
data = SpatialGaussian(data,[3 3],'kernel');

%reformat to the same size at motifs
temp = NaN(68,68,size(data,3));
temp(:,1:38,:)=data;
data = temp; clear temp;
data_flat = reshape(data,68*68,size(data,3));
data_flat(isnan(data_flat))=0;

%get only the pixels that are non masked in each
bad_pxl = (nanvar(data_flat,[],2)<=eps) + (nanvar(reshape(motifs,68*68,14*26),[],2)<=eps) >= 1;
data_flat(bad_pxl,:) = [];
motifs(bad_pxl,:,:) = [];

% H = H_orig([1,3,4,5,8,11],:);
H = H_orig; 
% labels = {'1','3','4','5','8','11'};
labels = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
% motifs = motifs(:,[1,3,4,5,8,11],:);
data = data_flat;
%Get the H onsets
Onset_times = cell(1,size(H,1));
Activity = H;
for j = 1:size(H,1)
   XHat = tensor_convolve(motifs(:,j,:),H(j,:));
   rho_frame = zeros(1,size(XHat,2));
   for i = 1:size(XHat,2)    
      rho_frame(i) = corr(XHat(:,i),data(:,i));
   end 
%     figure; hold on; plot(rho_frame); 
%     plot(H(j,:))
   [pks,locs] = findpeaks(rho_frame,'MinPeakHeight',0.4,'MinPeakDistance',1*15);
   Activity(j,:) = zeros(size(H(j,:)));
   Activity(j,locs) = pks;
   Onset_times{j} = locs;
end
set(gca,'xtick',(1:(10*60*15):size(H,2)))
set(gca,'xticklabel',floor(get(gca,'xtick')/15))
xlabel('Time (s)')
ylabel('Rho or Intensity')
fp.FormatAxes(gca)
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','ExampleCorrelation','Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all

%%
%compare with firing rates
temp  = load('firing_rate.mat','firing_rate_all');
fr = temp.firing_rate_all([1,2,4,5,6,7],:);
for i = 1:size(fr,1)
    fr(i,:) = smooth(fr(i,:),50);
end
fr = zscore(fr,0,2);%zscore
[coef, score, latent, ~, explained, mu] = pca(fr');
score = score';
cf = size(fr,2)/size(H,2);
data = cell(1,numel(Onset_times));
data_pc = cell(1,numel(Onset_times));
for i = 1:numel(Onset_times)
   temp = Onset_times{i};
   dur = floor(30*cf);
   snippet = NaN(size(fr,1),dur,numel(Onset_times));
   snippet_pc = NaN(size(score,1),dur,numel(Onset_times));
   for j = 1:numel(temp)
       try %avoid the end
          start = floor(temp(j)*cf-(15*cf));
          snippet(:,:,j) = fr(:,start:start+dur-1); 
          snippet_pc(:,:,j) = score(:,start:start+dur-1); 
       catch
       end
   end
   data{i} = snippet;
   data_pc{i} = snippet_pc;
end

%plot the average between the two
col = getColorPalet(size(fr,1));
figure('position',[28         558        1808         420]); hold on;
for i = 1:numel(data)
    subplot(1,numel(data),i); hold on; 
    x = 1:size(data{i},2);
    y = nanmean(data{i},3);
    err = nanstd(data{i},[],3)/sqrt(size(data{i},3));
    for j = 1:size(y,1)
        shadedErrorBar(x,y(j,:),err(j,:),'lineprops',{'color',col(j,:)});
    end
    title(sprintf('Motif %s',labels{i}),'FontWeight','normal');
    ylim([-0.7 1.3])    
    set(gca,'xticklabels',(-1:1:1))
    xlabel('time (s)')
    ylabel('Firing Rate (zscore)');
    fp.FormatAxes(gca);
end
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','FR','Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all

%PC Space
num_pcs = 3;
col = getColorPalet(size(fr,1));
figure('position',[28         558        1808         420]); hold on;
for i = 1:numel(data_pc)
    subplot(1,7,i); hold on; 
    x = 1:size(data_pc{i}(1:num_pcs,:,:),2);
    y = nanmean(data_pc{i}(1:num_pcs,:,:),3);
    err = nanstd(data_pc{i}(1:num_pcs,:,:),[],3)/sqrt(size(data_pc{i}(1:num_pcs,:,:),3));
    for j = 1:size(y,1)
        shadedErrorBar(x,y(j,:),err(j,:),'lineprops',{'color',col(j,:)});
    end
    title(sprintf('Motif %s',labels{i}),'FontWeight','normal');
    ylim([-0.7 1.3])    
    set(gca,'xticklabels',(-1:1:1))
    xlabel('time (s)')
    ylabel('Firing Rate (zscore)');
    fp.FormatAxes(gca);
end
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','PC_FR','Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all


%% Plot the average trajectory
fp = fig_params;
col = distinguishable_colors(14);
figure; hold on;
for i = 1:numel(data_pc)         
    pc_space = coef(:, 1:num_pcs)'*nanmean(data_pc{i}(:,:,:),3);    
      
%     pc_space_avg = movmean(pc_space_avg,10,2);

    a = linspace(75,255,size(pc_space,2));
    N = size(pc_space,2);      
    cdmap = uint8([repmat(col(i,:),size(pc_space,2),1)*255,a'])';
    
    if num_pcs ==3 
        p = plot3(pc_space(1,:),pc_space(2,:),pc_space(3,:),'linewidth',5);
    else
        p = plot(pc_space(1,:),pc_space(2,:),'linewidth',5);    
    end
    p.Color =[col(i,:),0.95]; drawnow;    
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cdmap); drawnow;
    
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
legend(labels{:});
grid on
fp.FormatAxes(gca);
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','PC Trajectory','Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all


% Simply plot an ANOVA comparing the average neuron activity at the center
cf = size(fr,2)/size(H,2);
data = cell(1,numel(Onset_times));
for i = 1:numel(Onset_times)
   temp = Onset_times{i};
   dur = floor(1*cf);
   snippet = NaN(size(fr,1),numel(Onset_times));
   for j = 1:numel(temp)
       try %avoid the end
          start = floor(temp(j)*cf-(0.5*cf));
          snippet(:,j) = nanmean(fr(:,start:start+dur-1),2); 
       catch
       end
   end
   data{i} = snippet;
end
data = cellfun(@(x) x(:,~isnan(sum(x))),data,'UniformOutput',0);
grp = cell(1,numel(data));
for i =1:numel(data)    
    grp{i} = ones(1,size(data{i},2))*str2num(labels{i});
end
grp = cat(2,grp{:});
temp = cat(2,data{:});

%Loop through the neurons 
figure('position',[28         558        1808         420]); hold on;
for i = 1:size(temp,1)
   p = anovan(temp(i,:)',{grp},'display','off');
   subplot(1,7,i); hold on; 
   boxplot(temp(i,:)',grp,'OutlierSize',0.1,'PlotStyle','compact')    
   title(sprintf('Neuron %d\np = %0.2g',i,p),'fontsize',14,'fontweight','normal');
   ylim([-2,4])
   fp.FormatAxes(gca);
end

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','FR_BoxPlotAnova','Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all

%% Now do the same thing for the power spectrums
%Get the one second of LFP for one channel after a given motif
for chan = [4,10,26,30]
%     chan = 4;
lfp = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\lfp.mat');
lfp = lfp.lfp(:,chan)';
cam = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\camera_signal.mat');
params = processing_params;
params.morlet_width = 5; %increase
params.mortlet_timepadding = 3;
params.morlet_fq = 2.^[1.5:0.1:6];
lfp_trim = AllignEphysToImaging(lfp,cam.camera_signal,[1,size(H,2)*2]);
lfp_wavelet = abs((MorletTransform(lfp_trim',params))').^2;


%not split into onset times
conversionfactor = size(lfp_trim,2)/size(H,2);
snippet = cell(max(cellfun(@(x) numel(x),Onset_times)),numel(Onset_times));
sample_freq = 1000;
L = 999;
for M =1:numel(Onset_times)
    M
    for i = 1:max(cellfun(@(x) numel(x),Onset_times))
        try %runs into issues if near end
            frames = floor(Onset_times{M}(i)*conversionfactor);
            snippet{i,M} = lfp_wavelet(:,frames-L:frames+L+1);

            lfp_temp = lfp_trim(1,frames-L:frames+L+1);
            N = length(lfp_temp);
            xdft = pmtm(lfp_temp); %pmtm used to be fft
            xdft = xdft(1:N/2+1);
            psdx = (1/(sample_freq*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            freq = 0:sample_freq/length(lfp_temp):sample_freq/2;
            pwr{i,M} = psdx';
        catch
            snippet{i,M} = NaN(46,L*2+2);
            pwr{i,M} = NaN(1,L+2);
        end
    end
end

avg_trace = {};
avg_freq = {};
for i = 1:numel(Onset_times)    
    avg_trace{i} = nanmean(cat(3,snippet{:,i}),3);
    avg_freq{i} = nanmean(cat(3,pwr{:,i}),3);
end
avg_trace = nanmean(cat(3,avg_trace{:}),3);
avg_freq = nanmean(cat(3,avg_freq{:}),3);


val = [-1, 1];
val2 = [0.33, 3];
for M = 1:numel(Onset_times)
    temp = nanmean(cat(3,snippet{:,M}),3).*params.morlet_fq';
%     temp = nanmean(cat(3,snippet{:,M}),3);
    
    %visualize average lfp for motif 1
    figure('position',[259 558 1278 420]) ; hold on;  
    subplot(1,2,1); axis square;
    sgtitle(sprintf('motif %s',labels{M}));
%     imagesc(log(temp),val); c=colorbar;
    imagesc(temp); c=colorbar;
    set(gca,'YDir','normal')
    ylabel(c,' power')
    ylim([0.5,numel(params.morlet_fq)+0.5]);
    xlim([0,size(temp,2)]);
    set(gca,'ytick',1:5:numel(params.morlet_fq),'yticklabel',round(params.morlet_fq(1:5:end),2))
    set(gca,'xtick',0:500:L*2+2,'xticklabel',0-(L+1):500:L+1);
    ylabel('frequency (Hz)');
    xlabel('time (ms)')
    axis square;

    colormap magma        
    %frequency spectrum
    temp = cellfun(@(x) smooth(x,5),pwr(:,M),'UniformOutput',0);
    temp_avg = nanmean(cat(3,temp{:}),3).*freq';
    temp_sem = nanstd(cat(3,temp{:}),[],3)/sqrt(numel(temp)).*freq';    
    subplot(1,2,2); axis square; hold on;    
    shadedErrorBar(freq,temp_avg,temp_sem,'lineprops',{'r','markerfacecolor','r'});
    plot(freq,temp_avg,'color',[0.75 0 0],'linewidth',2);
%     ylim(val2)
    xlim([0 75])
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    axis square;
    drawnow
 
end
fp = fig_params;

figure('position',[680   770   560   208]); hold on; 
plot(lfp_temp,'linewidth',2,'color','k')
xlabel('time (ms)');
ylabel('Amplitude (V)');
fp.FormatAxes(gca);
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg',sprintf('chan%doscillation_newnormalization',chan),'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Figures\LabMeeting7_2_2020',1); close all

end
























%% Plot the average fr during the peak
%compare with firing rates
temp  = load('firing_rate.mat','firing_rate_smooth');
fr = temp.firing_rate_smooth;
fr = zscore(fr,0,2); %zscore
[coef, score, latent, ~, explained, mu] = pca(fr');
score = score';
cf = size(fr,2)/size(H,2);
data = cell(1,numel(Onset_times));
data_pc = cell(1,numel(Onset_times));
for i = 1:numel(Onset_times)
   temp = Onset_times{i};
   dur = floor(1*cf);
   snippet = NaN(size(fr,1),numel(Onset_times));
   snippet_pc = NaN(size(score,1),numel(Onset_times));
   for j = 1:numel(temp)
       try %avoid the end
          start = floor(temp(j)*cf-(15*cf));
          snippet(:,j) = nanmean(fr(:,start:start+dur-1),2); 
          snippet_pc(:,j) = nanmean(score(:,start:start+dur-1),2); 
       catch
       end
   end
   data{i} = snippet;
   data_pc{i} = snippet_pc;
end

%remove empty 
num_pcs = 2;
data_pc = cellfun(@(x) x(:,~isnan(sum(x))),data_pc,'UniformOutput',0);
data = cellfun(@(x) x(:,~isnan(sum(x))),data,'UniformOutput',0);
grp = cell(1,numel(data));
for i =1:numel(data)    
    grp{i} = ones(1,size(data_pc{i},2))*i;
end
grp = cat(2,grp{:});
temp = cat(2,data{:});
[Y,loss] = tsne(temp');


pc = coef'*temp; 

gscatter(pc(1,:),pc(2,:),grp)




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%Get the one second of LFP for one channel after a given motif
lfp = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\lfp.mat');
lfp = lfp.lfp(:,4)';
cam = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\camera_signal.mat');
params = processing_params;
params.morlet_width = 5; %increase
params.mortlet_timepadding = 3;
params.morlet_fq = 2.^[1.5:0.1:6];
lfp_trim = AllignEphysToImaging(lfp,cam.camera_signal,[1,size(H,2)*2]);
lfp_wavelet = abs((MorletTransform(lfp_trim',params))').^2;

%not split into onset times
conversionfactor = size(lfp_trim,2)/size(H,2);
snippet = cell(max(cellfun(@(x) numel(x),Onset_times)),14);
sample_freq = 1000;
L = 1999;
X = 30;
for M =1:14
    M
    for i = 1:max(cellfun(@(x) numel(x),Onset_times))
        try %runs into issues if near end
            frames = floor([Onset_times{M}(i),Onset_times{M}(i)+X]*conversionfactor);
            snippet{i,M} = lfp_wavelet(:,frames(1):frames(1)+L);

            lfp_temp = lfp_trim(1,frames(1):frames(1)+L);
            N = length(lfp_temp);
            xdft = pmtm(lfp_temp); %pmtm used to be fft
            xdft = xdft(1:N/2+1);
            psdx = (1/(sample_freq*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            freq = 0:sample_freq/length(lfp_temp):sample_freq/2;
            pwr{i,M} = psdx';
        catch
            snippet{i,M} = NaN(46,L+1);
            pwr{i,M} = NaN(1,floor(L/2)+2);
        end
    end
end

%remove specific onsets so you can split by randonly bins
subset=0;
if subset ==1
    onset_mat = zeros(max(cellfun(@(x) numel(x),Onset_times)),14);
    for M = 1:14
        for i =  1:max(cellfun(@(x) numel(x),Onset_times))
            try
                onset_mat(i,M) = Onset_times{M}(i);
            end
        end
    end        
    times = [1,floor(size(H,2)/2)];
    idx = zeros(size(onset_mat));
    for M = 1:14
        temp = onset_mat(:,M);
        temp(temp<times(2) & temp>times(1))=1;
        idx(:,M) = temp;
    end
else
    idx = ones(max(cellfun(@(x) numel(x),Onset_times)),14);
end

% avg_trace = nanmean(cat(3,snippet{idx==1}),3);
% avg_freq = nanmean(cat(3,pwr{idx==1}),3);

avg_trace = {};
avg_freq = {};
for i = 1:14    
    avg_trace{i} = nanmean(cat(3,snippet{idx(:,i)==1,i}),3);
    avg_freq{i} = nanmean(cat(3,pwr{idx(:,i)==1,i}),3);
end
avg_trace = nanmean(cat(3,avg_trace{:}),3);
avg_freq = nanmean(cat(3,avg_freq{:}),3);


val = [0.5,2];
val2 = [0.5 2];
for M = 1:14
    temp = nanmean(cat(3,snippet{idx(:,M)==1,M}),3)./avg_trace;
    
    %visualize average lfp for motif 1
    figure('position',[259 558 1278 420]) ; hold on;  
    subplot(1,2,1); axis square;
    sgtitle(sprintf('motif %d',M));
    imagesc(temp,val); c=colorbar;
    set(gca,'YDir','normal')
    ylabel(c,'relative power')
    ylim([0.5,numel(params.morlet_fq)+0.5]);
    xlim([0.5,size(temp,2)]);
    set(gca,'ytick',1:5:numel(params.morlet_fq),'yticklabel',round(params.morlet_fq(1:5:end),2))
    set(gca,'xtick',0:500:L+1,'xticklabel',0:500:L+1);
    ylabel('frequency (Hz)');
    xlabel('time (ms)')
    axis square;

    colormap magma        
    %frequency spectrum
    temp = cellfun(@(x) smooth(x,5)./smooth(avg_freq,5),pwr(idx(:,M)==1,M),'UniformOutput',0);
    temp_avg = nanmean(cat(3,temp{:}),3);
    temp_sem = nanstd(cat(3,temp{:}),[],3)/sqrt(numel(temp));    
    subplot(1,2,2); axis square; hold on;    
    shadedErrorBar(freq,temp_avg,temp_sem,'lineprops',{'r','markerfacecolor','r'});
    plot(freq,temp_avg,'color',[0.75 0 0],'linewidth',2);
    ylim(val2)
    xlim([0 75])
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    axis square;
    drawnow
 
end
fp = fig_params;

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','chan4_shank1_oscillation_batch_all_2sec','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);close all


%% FIRING RATE DATA
%not split into onset times

reference_trace = load('Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse501_11_21_2019\Rec_Mouse501_112119\Block-1\reference_trace.mat');
cam = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\camera_signal.mat');
chan = [10,12,14,2,8,4,4];
shank = [1,1,2,2,1,2,1];
firing_rate_smooth = [];
firing_rate_all = [];
for i = 1:size(chan,2)%compile known good channels
    fprintf('\n working on %d',i);
    raw = TDTbin2mat('Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse501_11_21_2019\Rec_Mouse501_112119\Block-1\','STORE',['RAW',num2str(shank(i))],'CHANNEL',chan(i));
    raw = double(raw.streams.(['RAW',num2str(shank(i))]).data(1:numel(cam.camera_signal))); %written this way to trim extra 1 sample sometimes captured in ephys data vs time signal         
    %subtract reference
    raw = raw-reference_trace.reference_trace;
    ephys_trim = AllignEphysToImaging(raw,cam.camera_signal,[1,size(H,2)*2]);
    params = processing_params;
    params.s_spike_min_std = 10;
    [spike_wf, spike_times, spike_amp, ~, data_filt] = GetSpikes(ephys_trim,params);
    % figure; plot(data_filt(:,1:2.5e5)); hold on; plot(spike_times(1:100),-1*spike_amp(1:100),'marker','*','color','r','linestyle','none');

    %firing rate: number of spikes per 10 ms across whole recording
    edges = (1:params.sample_rate/100:size(ephys_trim,2));
    temp = discretize(spike_times,edges);
    temp(end)=[];
    fr = params.sample_rate ./diff(spike_times);
    time_bins = (1:1:numel(edges));
    firing_rate = time_bins; 
    for cur_bin = 1:numel(time_bins)
       firing_rate(time_bins(cur_bin)) = nanmean(fr(temp==time_bins(cur_bin)));
    end
    firing_rate(isnan(firing_rate))=0;
    firing_rate_smooth(i,:) = smooth(firing_rate,100);
    firing_rate_all(i,:) = firing_rate;
end
save('firing_rate.mat','firing_rate_smooth','firing_rate_all','chan','shank','-v7.3')


% firing_rate_smooth = movmean(firing_rate_smooth,100,2);
firing_rate_zscore = (firing_rate_smooth-nanmean(firing_rate_smooth,2))./std(firing_rate_smooth,[],2);
firing_rate_zscore = movmean(firing_rate_zscore,100,2);
conversionfactor = size(firing_rate_smooth,2)/size(H,2);
snippet = cell(max(cellfun(@(x) numel(x),Onset_times)),14);
sample_freq = 1000;
L = 1999;
X = 30;
for M =1:14
    for i = 1:max(cellfun(@(x) numel(x),Onset_times))
        try %runs into issues if near end
            frames = floor([Onset_times{M}(i),Onset_times{M}(i)+X]*conversionfactor);
            snippet{i,M} = firing_rate_zscore(:,frames(1):frames(1)+L);

        catch
            snippet{i,M} = NaN(size(chan,2),L+1);            
        end
    end
end

%remove specific onsets so you can split by randonly bins
subset=0;
if subset ==1
    onset_mat = zeros(max(cellfun(@(x) numel(x),Onset_times)),14);
    for M = 1:14
        for i =  1:max(cellfun(@(x) numel(x),Onset_times))
            try
                onset_mat(i,M) = Onset_times{M}(i);
            end
        end
    end        
    times = [1,floor(size(H,2)/2)];
    idx = zeros(size(onset_mat));
    for M = 1:14
        temp = onset_mat(:,M);
        temp(temp<times(2) & temp>times(1))=1;
        idx(:,M) = temp;
    end
else
    idx = ones(max(cellfun(@(x) numel(x),Onset_times)),14);
end

% zscore firing rate
[coef, score, latent, ~, explained, mu] = pca(firing_rate_zscore');
% [coef, score, latent, ~, explained, mu] = pca( squeeze(nanmean(cat(3,snippet{:}),2))');
% temp = squeeze(nanmean(cat(3,snippet{:}),2));
figure; hold on; 
col = distinguishable_colors(14);
label = {};
% get the average trace for each motif
COUNT = 1;
for M = [3,11]%,8,11]
    temp = cat(3,snippet{idx(:,M)==1,M});
%     temp = squeeze(nanmean(temp,2));
%     pc_space = coef(:, 1:3)'*temp;
%     
%     plot(pc_space(1,:),pc_space(2,:),'.','color',col(M,:),'markersize',2);
%     
    [x,y,z] = size(temp);
    temp = reshape(temp,[x,y*z]);
    pc_space = coef(:, 1:3)'*temp;  
    pc_space = reshape(pc_space,[size(pc_space,1),y,z]);
    pc_space_avg = nanmean(pc_space,3);    
      
    pc_space_avg = movmean(pc_space_avg,500,2);
%     plot the average trace in PC space for each motif    
    temp = pc_space_avg(:,1:1:size(pc_space_avg,2));
    a = linspace(75,255,size(temp,2));
    N = size(temp,2);      
    cd = uint8([repmat(col(COUNT,:),size(temp,2),1)*255,a'])';

    p = plot3(temp(1,:),temp(2,:),temp(3,:),'linewidth',5);
%     p = plot(temp(1,:),temp(2,:),'linewidth',5);    
    p.Color =[col(COUNT,:),0.95]; drawnow;    
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd); drawnow;
    label{COUNT} = {sprintf('Motif %d',M)};
    COUNT = COUNT+1;
    
%     plot(pc_space(1,:),pc_space(2,:),'.','color',col(COUNT,:),'markersize',1)
    
end
campos([-1.0505   -0.7537    0.3212]);
grid on; 
axis square
xlabel('PC1');ylabel('PC2');zlabel('PC3')

legend([label{:}]);
fp.FormatAxes(gca)
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-png','example_trajectory','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);
%% just the firing rate
% avg_trace = {};
% for i = 1:14    
%     avg_trace{i} = nanmean(cat(3,snippet{idx(:,i)==1,i}),3);
% end
% avg_trace = nanmean(cat(3,avg_trace{:}),3);

% % firing_rate_zscore = (firing_rate_smooth)./nanmean(firing_rate_smooth,2);
% conversionfactor = size(firing_rate_smooth,2)/size(H,2);
% snippet = cell(max(cellfun(@(x) numel(x),Onset_times)),14);
% sample_freq = 1000;
% L = 1999;
% X = 30;
% for M =1:14
%     for i = 1:max(cellfun(@(x) numel(x),Onset_times))
%         try %runs into issues if near end
%             frames = floor([Onset_times{M}(i),Onset_times{M}(i)+X]*conversionfactor);
%             snippet{i,M} = firing_rate_zscore(:,frames(1):frames(1)+L);
% 
%         catch
%             snippet{i,M} = NaN(size(chan,2),L+1);            
%         end
%     end
% end


val = [-0.5,0.5];
for M = [3,5,11]     
    temp = cellfun(@(x) x,snippet(idx(:,M)==1,M),'UniformOutput',0);
    temp_avg = nanmean(cat(3,temp{:}),3);
    
    %visualize average lfp for motif 1
    figure('position',[259 558 500 420]) ; hold on;      
    sgtitle(sprintf('motif %d',M));
    imagesc(temp_avg,val); c=colorbar;
    set(gca,'YDir','normal')
    ylabel(c,'relative power')
    ylim([0.5,size(temp_avg,1)+0.5]);
    xlim([0.5,size(temp_avg,2)]);
%     set(gca,'ytick',1:5:numel(params.morlet_fq),'yticklabel',round(params.morlet_fq(1:5:end),2))
    set(gca,'xtick',0:500:L+1,'xticklabel',0:500:L+1);
    ylabel('Firing Rate (normalized)');
    xlabel('time (ms)')
    axis square;
    colormap magma
end

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','firing_batch_all_2sec','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);close all

%% just the firing rate
avg_trace = {};
for i = 1:14    
    avg_trace{i} = nanmean(cat(3,snippet{idx(:,i)==1,i}),3);
end
avg_trace = nanmean(cat(3,avg_trace{:}),3);



% firing_rate_zscore = (firing_rate_smooth)./nanmean(firing_rate_smooth,2);
% fr_rate_single = movmean(firing_rate_smooth,100,2);
fr_rate_single = firing_rate_smooth;
conversionfactor = size(fr_rate_single,2)/size(H,2);
snippet = cell(max(cellfun(@(x) numel(x),Onset_times)),14);
sample_freq = 1000;
L = 1999;
X = 30;
for M =1:14
    for i = 1:max(cellfun(@(x) numel(x),Onset_times))
        try %runs into issues if near end
            frames = floor([Onset_times{M}(i),Onset_times{M}(i)+X]*conversionfactor);
            snippet{i,M} = fr_rate_single(:,frames(1):frames(1)+L);

        catch
            snippet{i,M} = NaN(size(chan,2),L+1);            
        end
    end
end

val = [-0.5,0.5];

figure('position',[259 558 500 420]) ; hold on; 
col = distinguishable_colors(14);
for M = [3,11]  
        
%     temp = cellfun(@(x) movmean(x,5,2),snippet(idx(:,M)==1,M),'UniformOutput',0);
    temp = cellfun(@(x) movmean(x,30,2),snippet(6,M),'UniformOutput',0);
    temp_avg = nanmean(cat(3,temp{:}),3);
    temp_sem = nanstd(cat(3,temp{:}),[],3)/sqrt(numel(temp));     
%     for i = 1:size(temp_avg,1)
% %         subplot(size(temp_avg,1),1,i)
%         shadedErrorBar(1:size(temp_avg,2),temp_avg(i,:),temp_sem(i,:));
%     end
    plot(1:size(temp_avg,2),temp_avg(2,:),'linewidth',2,'color',col(M,:));
%     ylim(val2)   
%     xlim([0.5,size(temp,2)]);
    set(gca,'xtick',0:500:L+1,'xticklabel',0:500:L+1);
    ylabel('Firing Rate (HZ)');
    xlabel('time (ms)')
%     axis square;    

end


handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','example_FR','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);close all
 
    
    
    %     %frequency spectrum
%     figure; hold on;
%     temp = cellfun(@(x) movmean(x,10,2),snippet(idx(:,M)==1,M),'UniformOutput',0);
%     temp_avg = nanmean(cat(3,temp{:}),3);
%     temp_sem = nanstd(cat(3,temp{:}),[],3)/sqrt(numel(temp));     
%     for i = 1:size(temp_avg,1)
%         shadedErrorBar(1:size(temp_avg,2),temp_avg(i,:),temp_sem(i,:));
%     end
% %     ylim(val2)   
% %     xlim([0.5,size(temp,2)]);
%     set(gca,'xtick',0:500:L+1,'xticklabel',0:500:L+1);
%     ylabel('Firing Rate (zscore)');
%     xlabel('time (ms)')
%     axis square; 
end
fp = fig_params;

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','oscillation_batch_all_2sec','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);close all
% % 
% % 
% % for M = [3,5,11,1,8]
% %     temp = cat(3,snippet{idx(:,M)==1,M});
% %     [x,y,z] = size(temp);
% %     temp = reshape(temp,[x,y*z]);
% %     pc_space = coef(:, 1:3)'*temp;  
% %     pc_space = reshape(pc_space,[size(pc_space,1),y,z]);
% %     pc_space_avg = nanmean(pc_space,3);    
% %       
% %     pc_space_avg = movmean(pc_space_avg,500,2);
% %     % plot the average trace in PC space for each motif    
% %     temp = pc_space_avg(:,1:1:size(pc_space_avg,2));
% %     a = linspace(25,255,size(temp,2));
% %     N = size(temp,2);
% % %     cd = [[linspace(1,col(M,1),N)',linspace(1,col(M,2),N)',linspace(1,col(M,2),N)']*255,ones(N,1)]';
% %        
% %     cd = uint8([repmat(col(M,:),size(temp,2),1)*255,a'])';
% % %     p = plot3(temp(1,:),temp(2,:),temp(3,:),'marker','none','linestyle','-','markersize',10);
% %     p = plot3(temp(1,:),temp(2,:),temp(3,:),'linewidth',8);
% %     p.Color(4) =0.75;
% %     set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
% % %     set(p.MarkerHandle,'FaceColorBinding','interpolated', 'FaceColorData',cd)
% % %     set(p.MarkerHandle,'EdgeColorBinding','interpolated', 'EdgeColorData',cd)
% %     
% %     
% %     cd=uint8([255,200,250,50,0; 0,50,250,150,200; 0,0,0,100,150; 179,150,200,70,50]);    
% %     h2b = plot(2:6,  15:-1:11, '.-r', 'LineWidth',8, 'DisplayName',' 0.7'); h2b.Color(4)=0.7;  % 30% transparen
% %     set(h2b.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
% % %     for i = 1:size(temp,2)
% % %        scatter3(temp(1,i),temp(2,i),temp(3,i),10,col(M,:),'o','filled','MarkerFaceAlpha',a(i),'MarkerEdgeAlpha',a(i)); 
% % %     end
% % end
% % set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
% % legend;

% % lfp_wave_all = {};
% % for M = 1:14
% % lfp_trim = cell(numel(Onset_times{M}),1);
% % lfp_wavelet = cell(numel(Onset_times{M}),1);
% % for i = 1:numel(Onset_times{M})
% %    frames = [Onset_times{M}(i)-30,Onset_times{M}(i)+30]*2;
% %    lfp_trim{i} = AllignEphysToImaging(lfp,cam.camera_signal,frames);
% %    lfp_wavelet{i} = abs((MorletTransform(lfp_trim{i}',params))').^2;  
% % end
% % lfp_trim = MakeCellsEqual(lfp_trim,2,0);
% % lfp_wavelet = MakeCellsEqual(lfp_wavelet,2,0);
% % lfp_wave_all{i,M} = lfp_wavelet;
% % 
% % %visualize average lfp for motif 1
% % temp = nanmean(cat(3,lfp_wavelet{:}),3).*params.morlet_fq';
% % figure; hold on;  title(sprintf('motif %d',M));
% % imagesc(temp); c=colorbar;
% % ylabel(c,'power (uV^2 * fq')
% % ylim([0.5,numel(params.morlet_fq)+0.5]);
% % xlim([0.5,size(temp,2)]);
% % set(gca,'ytick',1:5:numel(params.morlet_fq),'yticklabel',round(params.morlet_fq(1:5:end),2))
% % set(gca,'xtick',[0:500:4000],'xticklabel',[-2000:500:2000]);
% % ylabel('frequency (Hz)');
% % xlabel('time (ms)')
% % 
% % end


















