function [rectitle,ImgPath,ImgProbeLoc,EphysPath,motif_fits,FaceCam,BodyCam,behavioral_tp] = LoadDataDirectories_VR(cur_rec)
%Camden MacDowell = timeless
%returns a list of the imaging, ephys, behavioral, probe, etc. data
%directories for each recording from the 2021 recordings

if cur_rec == 1
    %331 - recording 1
    rectitle = 'Mouse 331 Recording 1';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_11_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_11_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_11_2021_RestingState_g1\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_11_2021\Cam_1_20210611-132100.avi';
    BodyCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_11_2021\Cam_0_20210611-132100.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});   
    motif_fits_header = '\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat';
    behavioral_tp = [573,340484]; %first frame with LED on and first with it OFF
elseif cur_rec == 2
    %331 - recording 2
    rectitle = 'Mouse 331 Recording 2';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse331_RestingState_NP_06_12_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_12_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse331_06_12_2021_RestingState_g1\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_12_2021\Cam_1_20210612-163959.avi';
    BodyCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_12_2021\Cam_0_20210612-163959.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse331_RestingState_NP_06_12_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});    
    motif_fits_header = '\w*Mouse331_RestingState_NP_06_11_2021\w*chunk\w*.mat';
    behavioral_tp=[510,340484];
elseif cur_rec == 3
    %332 - recording 1
    rectitle = 'Mouse 332 Recording 1';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_07_2021_RestingState_g1\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_1_20210607-144148.avi';
    BodyCam = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Cam_0_20210607-144148.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});    
    motif_fits_header = '\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat';
    behavioral_tp = [579,340478];
elseif cur_rec == 4
    %332 - recording 2
    rectitle = 'Mouse 332 Recording 2';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse332_RestingState_NP_06_08_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_08_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse332_06_08_2021_RestingState_g0\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_08_2021\Cam_1_20210608-124351.avi';
    BodyCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_08_2021\Cam_0_20210608-124351.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse332_RestingState_NP_06_08_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    motif_fits_header = '\w*Mouse332_RestingState_NP_06_07_2021\w*chunk\w*.mat';
    behavioral_tp = [577,340542]; %first frame with LED on and first with it OFF
elseif cur_rec == 5
    %334 - recording 1
    rectitle = 'Mouse 334 Recording 1';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse334_RestingState_NP_06_09_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_09_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_09_2021_RestingState_g0\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Cam_1_20210609-165130.avi';
    BodyCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Cam_0_20210609-165130.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});    
    motif_fits_header = '\w*Mouse334_RestingState_NP_06_09_2021\w*chunk\w*.mat';
    behavioral_tp = [575,340472];
elseif cur_rec == 6
    %334 - recording 2
    rectitle = 'Mouse 334 Recording 2';
    ImgPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging\Mouse334_RestingState_NP_06_10_2021_1dff_combined_processed.mat';
    ImgProbeLoc = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_10_2021_probe_coords.mat';
    EphysPath = '\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\ProcessedEphys\catgt_Mouse334_06_10_2021_RestingState_g1\ap_opts.mat';
    FaceCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_10_2021\Cam_1_20210610-122726.avi';
    BodyCam = '\\cup\buschman\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_10_2021\Cam_0_20210610-122726.avi';
%     [motif_fits,~] = GrabFiles('\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    motif_fits_header = '\w*Mouse334_RestingState_NP_06_10_2021\w*chunk\w*.mat';
    behavioral_tp=[572,340482];
else
    error('unknown rec number | limit is 6');     
end

%adapt to spock
if ~ispc
    ImgPath = ConvertToBucketPath(ImgPath);
    ImgProbeLoc = ConvertToBucketPath(ImgProbeLoc);
    EphysPath = ConvertToBucketPath(EphysPath);
    FaceCam = ConvertToBucketPath(FaceCam);
    BodyCam = ConvertToBucketPath(BodyCam);
    [motif_fits,~] = GrabFiles(motif_fits_header,0,{ConvertToBucketPath('\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits')});
    temp = motif_fits; %reorder so train-->test
    temp(1:2:end) = motif_fits(2:2:end);
    temp(2:2:end) = motif_fits(1:2:end);
    motif_fits = temp;
else
    [motif_fits,~] = GrabFiles(motif_fits_header,0,{'\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'});
    temp = motif_fits; %reorder so train-->test
    temp(1:2:end) = motif_fits(2:2:end);
    temp(2:2:end) = motif_fits(1:2:end);
    motif_fits = temp;
end


end %function