function Spock_CombineStacks(folder_path,save_dir)
    if ~ispc
        addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    end
    try
        [stack, opts] = CombineStacks(folder_path);
        [~,header] = fileparts(folder_path);
        fn = [save_dir header 'dff_combined'];
        save(fn,'stack','opts','-v7.3');
    catch
        gp = general_params; 
        warning('\nNo dffs with suffix %s were combined in %s',gp.dff_suffix, folder_path);
    end
end