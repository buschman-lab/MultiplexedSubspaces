classdef general_params  
    properties
        %Path Options
        local_bucket = 'Z:\';
        spock_bucket = '\jukebox\buschman\';
        repo_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\';
        dynamic_script_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\Spock\DynamicScripts\';
        
        dff_suffix = '.ome_dff.mat' %the suffix used to take individual dff and ID for combining e.g. '_dff.mat' for corrected recs and '_dff_uncorrected.mat' for uncorrected recs
        delete_singledff = 1; %flag to delete individual dffs after done concatenating
        
    end

    methods
 
    end

end











