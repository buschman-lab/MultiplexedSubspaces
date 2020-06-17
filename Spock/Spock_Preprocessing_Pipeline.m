function Spock_Preprocessing_Pipeline(in_fn,opts_path)
%% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

opts = load(opts_path);
opts = opts.prepro_log;

[path, fn] = fileparts(in_fn);

fprintf('\nPreprocessing... file %s',fn)
stack = PreProcess(in_fn,opts);

%save off data
fprintf('\n Saving data');
fn = [path filesep fn '_stack.mat'];
save(fn,'stack','opts','-v7.3');

fprintf('\n Complete');
    
end



% %DEPRECIATED AS OF 6/3/2020. NEED TO WAIT UNTIL YOU COMBINE TO DO THE HEMODYNAMIC CORRECTION TO MAKE
% %SURE INDICES ARE CORRECT!!!
% if numel(unique(opts.wavelength_pattern))>1 %if multiple wavelengths used
%    [dff, dff_b, ~] = HemodynamicCorrection(stack, opts); 
% else
%    fprintf('\n No hemodynamic correction');
%    dff_b = makeDFF(stack, opts); 
%    dff = [];
% end
% 
% fprintf('\n Saving data');
% %Save off corrected data if available
% if ~isempty(dff)
%     [path, fn] = fileparts(in_fn);
%     %save dff
%     fn = [path filesep fn '_dff.mat'];
%     save(fn,'dff','opts','-v7.3');
% end
%     
% %save off the uncorrected if desired
% if opts.save_uncorrected
%     [path, fn] = fileparts(in_fn);
%     dff = dff_b; %need to have it still named dff. 
%     fn = [path filesep fn '_dff_uncorrected.mat'];
%     save(fn,'dff','opts','-v7.3');
% end
