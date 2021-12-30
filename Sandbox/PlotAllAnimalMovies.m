savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RawDataVideos';
figtitle = 'Mouse 332 Recording 1';
ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined_processed.mat';
ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_probe_coords.mat';
EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
BehavPath = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_1_20210607-144148.avi';
% BehavPath2 = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_0_20210607-144148.avi';
preprocessedFLAG = 1; %if plotting the deconvolved or nondeconvolved data
behavstr = [579,340478];
% Rectangle on the behavioral video to crop to visualize eyeball
EyeCrop = [260,180,65,60];

CombinedEphysImagingVideo(ImgPath,BehavPath,EphysPath,ImgProbeLoc,EyeCrop,figtitle,savedir,behavstr,preprocessedFLAG)

try
ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined.mat';
CombinedEphysImagingVideo(ImgPath,BehavPath,EphysPath,ImgProbeLoc,EyeCrop,figtitle,savedir,behavstr,0)
catch
end

savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RawDataVideos';
figtitle = 'Mouse 334 Recording 1';
ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse334_RestingState_NP_06_09_2021_1dff_combined_processed.mat';
ImgProbeLoc = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_09_2021_probe_coords.mat';
EphysPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
BehavPath = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Cam_1_20210609-165130.avi';
BehavPath2 = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Cam_0_20210609-165130.avi';
preprocessedFLAG = 1; %if plotting the deconvolved (1) or nondeconvolved (anything else) data
behavstr = [575,340472]; %first frame on and last frame on + 1 
% Rectangle on the behavioral video to crop to visualize eyeball
EyeCrop = [260,195,65,60,0];%add fifth term if the video is preflipped (I fixed this in the code mid experiment so animal 332 is not flipped and 334 is flipped). 

try 
CombinedEphysImagingVideo(ImgPath,BehavPath,EphysPath,ImgProbeLoc,EyeCrop,figtitle,savedir,behavstr,preprocessedFLAG)
catch 
end
    
ImgPath = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_09_2021_1dff_combined.mat';
CombinedEphysImagingVideo(ImgPath,BehavPath,EphysPath,ImgProbeLoc,EyeCrop,figtitle,savedir,behavstr,0)


