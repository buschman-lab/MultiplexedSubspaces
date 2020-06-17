function Plot_ExampleEphysAndImaging()
%for RO1 grant with Tim

addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Ephys_Processing');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\FunctionsForLFP');
addpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Plotting');
addpath(genpath('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\GitHub_repository\Ephys_Imaging\Toolboxes\TDTSDK'));

frames = [4030,4045];
chan = 10;

%allign
lfp = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\lfp.mat');
cam = load('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Preprocessed\Rec_Mouse501_112119\camera_signal.mat');
lfp_trim = AllignEphysToImaging(lfp.lfp',cam.camera_signal,frames*2);


%preprocess imaging
gp = general_params;
gp.w_pca_denoise=0;
[data_norm,nanpxs] = ProcessAndSplitData('Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_Preprocessed\501_11_21_2019_1dff_combined.mat',[],gp);

%Framing of example
data = conditionDffMat(data_norm',nanpxs,[],[68,38,size(data_norm,2)]);

%get example ephys with spiking
reference_trace = load('Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse501_11_21_2019\Rec_Mouse501_112119\Block-1\reference_trace.mat');
reference_trace = reference_trace.reference_trace;
raw = TDTbin2mat('Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse501_11_21_2019\Rec_Mouse501_112119\Block-1','STORE',['RAW',num2str(1)],'CHANNEL',chan);
raw = double(raw.streams.(['RAW',num2str(1)]).data(1:numel(reference_trace))); %written this way to trim extra 1 sample sometimes captured in ephys data vs time signal         
%subtract reference
raw = raw-reference_trace;
params = processing_params;
[~, ~, ~,~,data_filt] = GetSpikes(raw,params);
data_filt_trim = AllignEphysToImaging(data_filt,cam.camera_signal,frames*2);

%plot the example imag
temp = data; 
temp(isnan(temp))=0;
data_smooth = imgaussfilt3(temp,[1 1 0.1]);
data_smooth = data_smooth(:,:,frames(1):frames(2));
for i = 1:size(data_smooth,3)
    figure; 
    imagesc(fliplr(data_smooth(:,:,i)),[0 0.15]);
    colormap magma
    axis equal; 
    axis off
    title(sprintf('%d',i));
end

%Ephys examples 
lfp = lfp_trim(chan,:)*100000;

%plot the filtered LFP
fp = fig_params;
figure; hold on; 
plot(lfp,'color','k');
set(gca,'xlim',[0,size(lfp,2)],'xtick',[0,size(lfp,2)]);
temp = get(gca,'xtick');
set(gca,'xtick',linspace(temp(1),temp(end),5),'xticklabel',linspace(0,1000,5));
ylabel('Amplitude (uV)');
xlabel('time (ms)');
set(gca,'position',[0.1300    0.2850    0.7750    0.6400]);
set(gcf,'position',[ 801   382   685   213]);
fp.FormatAxes(gca);


figure; hold on; 
plot(data_filt_trim*100000,'color','k');
set(gca,'xlim',[0,size(data_filt_trim,2)],'xtick',[0,size(data_filt_trim,2)]);
temp = get(gca,'xtick');
set(gca,'xtick',linspace(temp(1),temp(end),5),'xticklabel',linspace(0,1000,5));
ylabel('Amplitude (uV)');
set(gca,'ylim',[-60 20]);
xlabel('time (ms)');
set(gca,'position',[0.1300    0.2850    0.7750    0.6400]);
set(gcf,'position',[ 801   382   685   213]);
fp.FormatAxes(gca);

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','example_ephysimaging_smooth','Z:\Users\Camden\GrantsWithTim\Ro1_2020',1);






















end %function








