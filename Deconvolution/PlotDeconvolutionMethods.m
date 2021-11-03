function PlotDeconvolutionMethods(file_path,save_path,dff,varargin)
%Camden MacDowell - timeless
%file_path goes to a preprocessed (i.e. dff) file with dff, probe_coords, opts, and nanpxs
%time window is frames to process and plot
opts.time_window = [800,2868]; %right now only works if both even or both odd
opts.video_flag = 1; 
opts.timebin = 1; 
opts.smoothkernel = [3 3]; %size in pixels
opts.pixel_scale = 2; %1=scale each pixel trace between zero and one. 2=zscore. 3=nothing

opts = ParseOptionalInputs(opts,varargin);

% file_path = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined.mat';
gp = general_params_corticaldynamics;

%load data
if isempty(dff)
    dff = load(file_path);
end
data = dff.dff(opts.time_window(1):opts.time_window(2),:);

%GLM | load trained GLM (see GenerateGLM_spock.sh)
load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging' filesep gp.fit_glm],'trained_opts');
temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list');
mdl_indx = find(strcmp(temp.dff_list,file_path)==1); %get the neural network that corresponds to this animal
kernel = flipud(trained_opts{mdl_indx}.glmkernel); %#ok<*FNDSB>
%loop through each pixel
data_glm = data;
for px = 1:size(data,2)            
    data_glm(:,px) = convn(padarray(data(:,px)',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid');    
end

%lucy richardson
data_lr = data;
for px = 1:size(data,2)
   data_lr(:,px) = lucric(data(:,px),0.89,1,size(kernel,1)-1);
end

%fNN | load trained fNN (see GeneratefNN_spock.sh)        
load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list');
net_indx = find(strcmp(temp.dff_list,file_path)==1); %get the neural network that corresponds to this animal
params = trained_opts{net_indx}.feedforwardparams;
net =trained_opts{net_indx}.shallowfeedforward;
%loop through each pixel
data_fnn = data;
for px = 1:size(data,2)            
    xtemp = createRollingWindow(data(:,px), params.win)'; %t-n:t-1        
    stPred = net(xtemp)';        
    %NaN pad to match timepoints
    data_fnn(:,px)=NaN;
    data_fnn(ceil(params.win/2):end-floor(params.win/2),px)=stPred;
end


%remove pad in all data
padidx = sum(isnan(data_fnn),2)==size(data_fnn,2); 
data_fnn(padidx,:)=[];
data_glm(padidx,:)=[];
data_lr(padidx,:)=[];
data(padidx,:)=[];


%% optional post processing

%scale within each pixel (to compensate for difference from the craniotomy)
if opts.pixel_scale==1
    for i = 1:size(data_fnn,2)
        data_fnn(:,i) = (data_fnn(:,i)-min(data_fnn(:,i)))/(prctile(data_fnn(:,i),99)-min(data_fnn(:,i)));
        data_glm(:,i) = (data_glm(:,i)-min(data_glm(:,i)))/(prctile(data_glm(:,i),99)-min(data_glm(:,i)));
        data(:,i) = (data(:,i)-min(data(:,i)))/(prctile(data(:,i),99)-min(data(:,i)));
        data_lr(:,i) = (data_lr(:,i)-min(data_lr(:,i)))/(prctile(data_lr(:,i),99)-min(data_lr(:,i)));
    end
elseif opts.pixel_scale==2 %zscore
    for i = 1:size(data_fnn,2)
        data_fnn(:,i) = zscore(data_fnn(:,i));
        data_glm(:,i) = zscore(data_glm(:,i));
        data(:,i) = zscore(data(:,i));
        data_lr(:,i) = zscore(data_lr(:,i));
    end    
end

if opts.timebin %bin in time
    data_fnn = data_fnn(1:2:end,:)+data_fnn(2:2:end,:);
    data_glm = data_glm(1:2:end,:)+data_glm(2:2:end,:);
    data_lr = data_lr(1:2:end,:)+data_lr(2:2:end,:);
    data = data(1:2:end,:)+data(2:2:end,:);
end

%condition
data_fnn = conditionDffMat(data_fnn,dff.nanpxs);
data_glm = conditionDffMat(data_glm,dff.nanpxs);
data_lr = conditionDffMat(data_lr,dff.nanpxs);
data = conditionDffMat(data,dff.nanpxs);

[x,y,z] = size(data);
mask = isnan(reshape(data(:,:,1),[x^2,1]));

if ~isempty(opts.smoothkernel) %spatially smoothing for visualization
    data(isnan(data))=0;
    data_fnn(isnan(data_fnn))=0;
    data_glm(isnan(data_glm))=0;
    data_lr(isnan(data_lr))=0;
    
    for n = 1:size(data,3)
        data(:,:,n) = imgaussfilt(data(:,:,n),'filtersize',opts.smoothkernel);
        data_fnn(:,:,n) = imgaussfilt(data_fnn(:,:,n),'filtersize',opts.smoothkernel);
        data_glm(:,:,n) = imgaussfilt(data_glm(:,:,n),'filtersize',opts.smoothkernel);
        data_lr(:,:,n) = imgaussfilt(data_lr(:,:,n),'filtersize',opts.smoothkernel);        
    end
    
    %remask to NaN
    data = reshape(data,[x^2,z]); data(mask,:) = NaN; data = reshape(data,[x,y,z]);
    data_fnn = reshape(data_fnn,[x^2,z]); data_fnn(mask,:) = NaN; data_fnn = reshape(data_fnn,[x,y,z]);
    data_glm = reshape(data_glm,[x^2,z]); data_glm(mask,:) = NaN; data_glm = reshape(data_glm,[x,y,z]);
    data_lr = reshape(data_lr,[x^2,z]); data_lr(mask,:) = NaN; data_lr = reshape(data_lr,[x,y,z]);    
end


%% plotting
%first plot one frame with the probe locations/craniotomy ring
%plot an example distribution for each method 
clear Frames
if opts.pixel_scale==1
    close all ; figure('position',[14 559 1864 420]);
    for i = 1:size(data,3)
        cla
        subplot(1,4,1); imagesc(fliplr(data(:,:,i)),[0 2]); colormap magma; title('original'); axis off; axis square; colorbar
        subplot(1,4,2); imagesc(fliplr(data_lr(:,:,i)),[0 2]); colormap magma; title('LR'); axis off; axis square; colorbar
        subplot(1,4,3); imagesc(fliplr(data_glm(:,:,i)),[0 2]); colormap magma; title('GLM'); axis off; axis square; colorbar
        subplot(1,4,4); imagesc(fliplr(data_fnn(:,:,i)),[0 2]); colormap magma; title('fNN'); axis off; axis square; colorbar     
        sgtitle(sprintf('frame %d',i));
        pause(0.01);
        Frames(i) = getframe(gcf);
    end
elseif opts.pixel_scale==2
    close all ; figure('position',[14 559 1864 420]);
    for i = 1:size(data,3)
        cla
        subplot(1,4,1); imagesc(fliplr(data(:,:,i)),[-3 3]); colormap magma; title('original'); axis off; axis square; colorbar
        subplot(1,4,2); imagesc(fliplr(data_lr(:,:,i)),[-3 3]); colormap magma; title('LR'); axis off; axis square; colorbar
        subplot(1,4,3); imagesc(fliplr(data_glm(:,:,i)),[-3 3]); colormap magma; title('GLM'); axis off; axis square; colorbar
        subplot(1,4,4); imagesc(fliplr(data_fnn(:,:,i)),[-3 3]); colormap magma; title('fNN'); axis off; axis square; colorbar     
        sgtitle(sprintf('frame %d',i));
        pause(0.01);
        Frames(i) = getframe(gcf);
    end    
else
    close all ; figure('position',[14 559 1864 420]);
    for i = 1:size(data,3)
        cla
        subplot(1,4,1); imagesc(fliplr(data(:,:,i)),round([prctile(data(:),1),prctile(data(:),99)])); colormap magma; title('original'); axis off; axis square; colorbar
        subplot(1,4,2); imagesc(fliplr(data_lr(:,:,i)),round([prctile(data_lr(:),1),prctile(data_lr(:),99)])); colormap magma; title('LR'); axis off; axis square; colorbar
        subplot(1,4,3); imagesc(fliplr(data_glm(:,:,i)),round([prctile(data_glm(:),1),prctile(data_glm(:),99)])); colormap magma; title('GLM'); axis off; axis square; colorbar
        subplot(1,4,4); imagesc(fliplr(data_fnn(:,:,i)),round([prctile(data_fnn(:),1),prctile(data_fnn(:),99)])); colormap magma; title('fNN'); axis off; axis square; colorbar     
        sgtitle(sprintf('frame %d',i));
        pause(0.01);
        Frames(i) = getframe(gcf);
    end
end
close; 
% create the video writer 
writerObj = VideoWriter([save_path '.avi']);
writerObj.FrameRate = 5;
open(writerObj);
for i=1:length(Frames)
    writeVideo(writerObj, Frames(i));
end
close(writerObj);


%% video 




