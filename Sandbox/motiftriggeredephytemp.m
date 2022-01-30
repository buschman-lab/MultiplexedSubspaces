%Goal is ephys split into trials according to motif occurance

%data
ephys_list = dir('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys'); 
ephys_list(1:2)=[]; %remove '.' and '..'
ephys_list(end)=[];
ephys_list=arrayfun(@(n) fullfile(ephys_list(n).folder,ephys_list(n).name),1:size(ephys_list,1),'UniformOutput',0);
[ap_fn,~] = GrabFiles('ap_opts.mat',0,ephys_list);
[motif_fits,~] = GrabFiles('\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
[orig_data,~] = GrabFiles('\w*processed\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging'});

%%todo
%order the motif onset chunks (probable test train)


%for each recording, epoch, get the onsets and concatenate 

%output is going to be recording x probe cell aray containing neuron,
%motif, time tensor of spike times. 



%time, motif, probe, neuron, recording
onsets = cat(1,onsets{:});
%add the duration of each chunk
temp = onsets;
onsets(1:2:88,:) = temp(2:2:88,:);
onsets(2:2:end,:) = temp(1:2:end,:);

for i = 2:size(onsets,1)
    temp = onsets(i,:);
    temp = cellfun(@(x) x+(900*(i-1)),temp,'UniformOutput',0);
    onsets(i,:) = temp;
end
motif_onset = [onsets{:,1}];

%plot the mean activity in the four craniotomies in the raw data
%plot the mean firing rate for the cortical region

%load ephys
p1 = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\AP_Probe1.mat');
p2 = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\AP_Probe2.mat');
p3 = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\AP_Probe3.mat');
p4 = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\AP_Probe4.mat');

spike_opts_list = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
%load imaging
dff = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_11_2021_1dff_combined_processed.mat','nanpxs','data_norm','gp','opts');
%need to move and change the probes a bit
coords = cellfun(@(x) x(1,:), probe_coords, 'UniformOutput', false);

%bin ephys
temporal_bin = 1;

%shift onsets to account for deconvolution data loss 
[st_mat_super,~] = LoadSpikes(spike_opts_fn,'depth',[0 600]);
[st_mat,~,st_depth] = LoadSpikes(spike_opts_fn);
[trace,~] = dffROI(dff,nanpxs,coords);

%get the imaging traces 
dff_roi = ROItrace(dff);

fp = fig_params;

w_good = W_basis;
w_good(:,[3,10],:) = [];
w_good = MaskTensor(w_good,nanpxs,1);

%% motif triggered figures
% plot the cortical area activity
for  n=1:8  
close all
motif_onset = [onsets{:,n}];
temp = arrayfun(@(x) dff_roi(x-5:x+35,:),motif_onset,'UniformOutput',0);
temp = cat(3,temp{:});
col = getColorPalet(4);
figure; hold on; 
for i =1:4
   shadedErrorBar(1:size(temp,1),nanmean(temp(:,i,:),3),sem(temp(:,i,:),3),'lineprops',{'color',col(i,:)});
end
ylabel({'mean zscore','in craniotomy'})
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
title(sprintf('deconvolved dff, motif %d',n));

% motif activity
[trace,~] = dffROI(squeeze(w_good(:,n,:))',nanpxs,coords);
figure; hold on; 
for i = 1:4
    plot(cat(1,zeros(5,1),trace(:,i)),'color',col(i,:),'linewidth',2);
end
ylabel({'motif activity','in craniotomy'})
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
title(sprintf('Motif dff, motif %d',n));

% plot the average ephys activity in superficial probe
temp = cellfun(@(x) nanmean(x,2), st_mat_super,'UniformOutput',0);
temp = [temp{:}];
temp = temp./std(temp,[],1);
temp = arrayfun(@(x) temp(x-5:x+35,:),motif_onset,'UniformOutput',0);
temp = cat(3,temp{:});
figure; hold on; 
for i =1:4
   shadedErrorBar(1:size(temp,1),nanmean(temp(:,i,:),3),sem(temp(:,i,:),3),'lineprops',{'color',col(i,:)});
end
ylabel({'std-normalized activity'})
title(sprintf('top 600mm, mean spiking, motif %d',n));
line([5,5],get(gca,'ylim'),'linewidth',2,'color','r');
% plot the motif heatmap organized by depth
temp_mat = cellfun(@(x) zscore(x), st_mat,'UniformOutput',0);
%reorganize by dpeth 
for i = 1:numel(temp_mat)
   [~,idx] = sort(st_depth{i},'ascend');
   temp_mat{i} = temp_mat{i}(:,idx);   
end

figure('position',[108 391 1654 588]); hold on; 
label = {'par','ss','v-thal','m2-pfc'};
for i = 1:numel(temp_mat)
    subplot(1,4,i);
    temp = arrayfun(@(x) temp_mat{i}(x-5:x+35,:),motif_onset,'UniformOutput',0);
    temp = nanmean(cat(3,temp{:}),3);
    imagesc(temp'); colormap gray; colorbar
    title(sprintf('%s \ndepth=%0.2f cm',label{i},max(st_depth{i})/1000))
    sgtitle(sprintf('motif %d',n));
    line([5,5],get(gca,'ylim'),'linewidth',2,'color',col(i,:));
    ylabel('neuron');
end

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\motif_trig_fig';
saveCurFigs(get(groot, 'Children'),'-dpng',sprintf('TriggeredActivity_motif%d',n),savedir,0); close all
end



