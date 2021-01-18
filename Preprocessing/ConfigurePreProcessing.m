function opts = ConfigurePreProcessing(varargin)
opts.fixed_image = 'first'; %image to use as reference. first, middle, mean, variance, max
opts.crop_h = 540; %height of croped image in pixels. 
opts.crop_w = 540; %width of cropped image in pixels.
opts.x_bregma_margin = 280;  %For outlininign brain: Number of pixels to anterior to bregma (default = 280)
opts.y_bregma_margin = 230; %for outlining brain: number of pixels to keep lateral to bregma (default - 230)
opts.mask_vasc = 1; %Mask vasculature
opts.manual_mask = 1; %Add additional manual masking to image (e.g. if imperfection in skull); 
opts.vasc_std = 1.5; %number of standard devidations from mean of ref_img to consider vasculature (2.5 works well)
opts.close_disk_size = 2; %pixel size of disk used to remove salt/pepper noise on vascmask. 
opts.verbose = 1; %how chatty do we want to be
opts.mask_brain_outline = 1; %Create mask for outside of brain (manually made). 
opts.mask_brain_outline_dir = [fileparts(which('ConfigurePreProcessing.m')) filesep 'brainoutline.mat']; 
opts.save_uncorrected = 0; %save the uncorrected stacks? 
opts.spatial_bin_factor = 8;
opts.method = 'movingavg'; %'movingavg','mean','median','mode'
opts.method_window = 30; %moving window duration (in s)
opts.fps = 15; %the frame rate of each wavelength (integer)
opts.wavelength_pattern = [1,2]; %excitation wavelength sequence used
opts.correction_wavelength = 2; %wavelength to be used for correction
opts.detrend = 0; %also include linear detrending(not needed if DFF)

opts = ParseOptionalInputs(opts,varargin);
end