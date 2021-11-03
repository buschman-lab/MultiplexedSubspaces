function Spock_Retinotopy_Pipeline(in_fn,opts_path)
%% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
end

opts = load(opts_path);
opts = opts.prepro_log;

[path, fn] = fileparts(in_fn);

%get trial index from the filename in case
fn=erase(fn,'_MMStack_Pos0.ome');
index = regexp(fn, '_', 'split');
index = str2num(index{end});

fprintf('\nPreprocessing... file %s',fn)
stack = PreProcess(in_fn,opts);

%use the same duration baseline period for all recordings
baselinedur = min(opts.frame_index(:,1)); 
if mod(baselinedur,2)==1 %if odd, then subtract 1
    baselinedur=baselinedur-1;
end

%change opts.window to the baseline period
opts.method_window = [opts.frame_index(index,1)-baselinedur+1,opts.frame_index(index,1)-1]; %substract 1 because this is the first frame AFTER a stim onset

if numel(unique(opts.wavelength_pattern))>1 %if multiple wavelengths used
   [dff, dff_b, ~] = HemodynamicCorrection_Sensory(stack, opts); 
   %ImpactOfHemoCorrection(dff,dff_b,dff_h)
else
   fprintf('\n No hemodynamic correction');
   dff_b = makeDFF(stack, opts); 
   dff = [];
end

%trim to the stimulus period
dff = dff(:,:,opts.frame_index_blue(index,1)-1:opts.frame_index_blue(index,2));
dff_b = dff_b(:,:,opts.frame_index_blue(index,1)-1:opts.frame_index_blue(index,2));

%save off data
fprintf('\n Saving data');
fn = [path filesep fn '_stack.mat'];
save(fn,'dff','dff_b','opts','index','-v7.3');

fprintf('\n Complete');
    
end


