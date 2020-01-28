%QUALITY CHECK FOR IMAGING DATA
%configure preprocessing options
opts.fixed_image = 0; %image to use as reference. 0=first, 1, middle, 2=mean
opts.crop_h = 540; %height of croped image in pixels. 
opts.crop_w = 540; %width of cropped image in pixels.
opts.x_bregma_margin = 280;  %For outlininign brain: Number of pixels to anterior to bregma (default = 280)
opts.y_bregma_margin = 230; %for outlining brain: number of pixels to keep lateral to bregma (default - 230)
opts.mask_vasc = 0; %Mask vasculature
opts.manual_mask = 0; %Add additional manual masking to image (e.g. if imperfection in skull); 
opts.vasc_std = 1.5; %number of standard devidations from mean of ref_img to consider vasculature (2.5 works well)
opts.close_disk_size = 2; %pixel size of disk used to remove salt/pepper noise on vascmask. 
opts.verbose = 1; %how chatty do we want to be
opts.mask_brain_outline = 1; %Create mask for outside of brain (manually made). 
opts.mask_brain_outline_dir = [fileparts(which('ConfigurePreProcessing.m')) filesep 'brainoutline_hemi.mat']; 
opts.save_uncorrected = 1; %save the uncorrected stacks? 
opts.spatial_bin_factor = 4;
opts.method = 'mean'; %'movingavg','mean','median','mode'
opts.fps = 15; %the frame rate of each wavelength (integer)
opts.wavelength_pattern = [1,2]; %excitation wavelength sequence used
opts.correction_wavelength = 2; %wavelength to be used for correction
opts.detrend = 1; %also include linear detrending 

%select folders to process and grab the first file from each rec.
% file_list = 'INSERT FILE PATH HERE';

%EXAMPLE OF A GOOD RECORDING: 'Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse501_11_21_2019\501-11-21-2019_1\501-11-21-2019_1_MMStack_Pos0.ome.tif';
%EXAMPLE OF A RECORDING WHERE NECK HAIR GOT IN WAY:'Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse495_1_26_2020\Mouse495_1_26_2020\495_1_26_2020_1\495_1_26_2020_1_MMStack_Pos0.ome.tif';
file_name = 'Z:\Rodent Data\Wide Field Microscopy\ParietalCortex_Ephys_Widefield\Mouse495_1_27_2020\495_1_27_2020_1\495_1_27_2020_1_MMStack_Pos0.ome.tif';

%Grab reference image  
ref_imgs = GetReferenceImage(file_name,opts.fixed_image);

%manual allignment 
prepro_log = ManualAlignment(ref_imgs,opts);
prepro_log.tform = [];
prepro_log.mask = [];

%make dff 
warning ('off','all'); %Suppress missing tiff metadata warning
img_count = 500; %just grab first 500 images

%preallocate stack
stack = NaN(floor(opts.crop_w/opts.spatial_bin_factor),...
    floor(opts.crop_h/opts.spatial_bin_factor),...
    img_count);

%read tiff
in_tiff = Tiff(file_name, 'r');

%Loop through each image
for cur_img_ind = 1:img_count
    in_tiff.setDirectory(cur_img_ind);   
    cur_img = in_tiff.read;
    
    %Preprocess image
    cur_img = PreProcessImage(cur_img,prepro_log);

    %Add to stack
    stack(:,:,cur_img_ind) = cur_img; 
    
    %Chatty 
    if prepro_log.verbose
        if mod(cur_img_ind,round(0.02*img_count)) ==0
            fprintf('\t%g%% Complete\n', round(cur_img_ind./img_count*100,2));
        end
    end
    
    
end %image loop
warning ('on','all');

%%
[dff, dff_b, ~] = HemodynamicCorrection(stack, prepro_log);
dff_b(isnan(dff_b))=0;
dff(isnan(dff))=0;

% Filter dff stack
Fs = 15;  % Sampling Frequency
N     = 11;   % Order
BW    = 0.2;  % Bandwidth
Apass = 5;    % Bandwidth Attenuation
filtband1 = [0.1*2/Fs 4*2/Fs];
[B,A] = butter(10,[filtband1(1),filtband1(2)],'bandpass');
[nX,nY,nZ] = size(dff_b);
filtStack_b = zeros(nX,nY,nZ);
filtStack = zeros(nX,nY,nZ);
for xx = 1:nX
    for yy = 1:nY
      filtStack_b(xx,yy,:) = filtfilt(B,A,squeeze(dff_b(xx,yy,:)))+nanmean(dff_b(xx,yy,:));% add back dc
      filtStack(xx,yy,:) = filtfilt(B,A,squeeze(dff(xx,yy,:)))+nanmean(dff(xx,yy,:));% add back dc
    end
    if mod(xx,round(0.1*nX)) ==0
       fprintf('\t%g%% done filtering...\n', xx./nX*100);
    end
end
dff_b = filtStack_b;
dff = filtStack;
%%
figure('position',[51         309        1651         669]);
mind = 0;
maxd = 2;
colormap magma
idx = 0;
s1 = subplot(1,2,1); 
s2 = subplot(1,2,2); 
for i = 1:size(dff_b,3)
    axes(s1); axis square; axis off; set(gca,'ydir','reverse');
    imagesc(dff_b(:,:,i));
    caxis([mind maxd]); 
    
    axes(s2); axis square; axis off; set(gca,'ydir','reverse');
    imagesc(dff(:,:,i));
    axis off
    caxis([mind maxd]); 
    sgtitle(sprintf('frame %d,',i))
    drawnow
    pause(0.01)
    hold off
end


