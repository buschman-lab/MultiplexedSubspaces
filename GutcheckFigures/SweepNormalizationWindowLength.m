function SweepNormalizationWindowLength(mouse)
%Camden MacDowell - timeless
%Cycles through recordings and sweeps the normalization window length.
%Saves off for subsequent use on
%'CompareDeconvolution_NormalizationWindowLength'
% mouse = {'331_1','331_2','332_1','332_2','334_1','334_2'}
fp = fig_params_deconvolutionpaper;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Normalization\';
if ~exist(savedir,'dir'); mkdir(savedir); end
switch mouse
    case '332_1'
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1'}); %cranio raw stack
    case '332_2'        
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse332_06_08_2021\Mouse332_RestingState_NP_06_08_2021_1'}); %cranio raw stac
    case '331_1'
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_11_2021\Mouse331_RestingState_NP_06_11_2021_1'}); %cranio raw stac        
    case '331_2'
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse331_06_12_2021\Mouse331_RestingState_NP_06_12_2021_1'}); %cranio raw stac        
    case '334_1'
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_09_2021\Mouse334_RestingState_NP_06_09_2021_1'}); %cranio raw stac  
    case '334_2'        
        craniotomy_file = GrabFiles('_Pos0.ome_stack.mat',0,...
            {'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\RestingState_Neuropixels\Mouse334_06_10_2021\Mouse334_RestingState_NP_06_10_2021_1'}); %cranio raw stac 
end

processed_path = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\'; %for probe coords
[~,fn_header] = fileparts(craniotomy_file{1});
probe_coords = load([processed_path,erase(fn_header,'_MMStack_Pos0.ome_stack'),'dff_combined.mat'],'probe_coords'); 
probe_coords=probe_coords.probe_coords;
header = erase(fn_header,'_MMStack_Pos0.ome_stack');

%% load the data and save off
[stack,opts] = GrabAFewMinutes(craniotomy_file,0); %dff_flag==1-->hemo 2=nohemo 0=nothing

opts.method = 'zscore';

% sweep method window (in seconds)
w = [5,10,30,60,120,300,600,1200];

for i = 1:numel(w)
    fprintf('\n\t Working on mouse %s window %d',mouse,w(i));
    opts.method_window = w(i);
    dff = makeDFF(stack, opts, 'dff', opts.method_window);
    [dff,nanpxs] = conditionDffMat(dff);
    save([savedir header sprintf('_window%d_seconds.mat',w(i))],'opts','w','probe_coords','dff','nanpxs','-v7.3');
end


end %main function 

function [stack,opts] = GrabAFewMinutes(file_list,dff_flag)
    %grab multiple preprocessed stacks (but not the full hour to save time)
    for i = 1:39 %five should be enough its  a little over 20 min at 30Hz... there should be ~40 files in total for the 90 minute recs
       [file_path,fn] = fileparts(file_list{1}); 
       if i==1
           temp_stack = load(file_list{1},'stack','opts');
           opts = temp_stack.opts;
           stack = temp_stack.stack;
       else
           temp_stack = load([file_path filesep erase(fn,'.ome_stack') sprintf('_%d.ome_stack.mat',i-1)],'stack');
           stack = cat(3,stack,temp_stack.stack);
       end
    end
    
    if dff_flag==1 %corrected
       stack = HemodynamicCorrection(stack, opts);
    elseif dff_flag==2 %uncorrected
        [~,stack] = HemodynamicCorrection(stack, opts);       
    else %filter to remove
        stack = filterstack(stack, 30, [0.1 4], 'lowpass', 1, 0);
    end
    
    mask = repmat(imresize(opts.mask,[68 68]),1,1,size(stack,3)); 
    stack(mask==0)=0; %correct in case this is yet to be included in the preprocessing
end