function Spock_CombineStacks(folder_path,save_fn)
    if ~ispc
        addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    end
    try
        [stack, opts] = CombineStacks(folder_path);
        save(save_fn,'stack','opts','-v7.3');
        
%         %make figures
%         fh = DetectDictalEvents(save_fn);
    catch
        gp = general_params; 
        warning('\nNo dffs with suffix %s were combined in %s',gp.dff_suffix, folder_path);
    end
end

