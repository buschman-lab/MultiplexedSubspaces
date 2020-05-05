classdef general_params  
    properties
        %Path Options
        local_bucket = 'Z:\';
        spock_bucket = '\jukebox\buschman\';
        repo_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\';        
        dynamic_script_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\Spock\DynamicScripts\';
        processing_intermediates = 'Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Scratch_Processing_Intermediates\'; %location of the intermediate files in the processing pipeline
        
        %dynamic script options
        sbatch_time = 59;
        sbatch_exclude = 'redshirt-n[12-49]';
        sbatch_memory = 16;
        sbatch_matlabversion = 'R2018a';
        sbatch_path = "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/";
        sbatch_name = [];
        
        
        dff_suffix = '.ome_dff.mat' %the suffix used to take individual dff and ID for combining e.g. '_dff.mat' for corrected recs and '_dff_uncorrected.mat' for uncorrected recs
        delete_singledff = 1; %flag to delete individual dffs after done concatenating
        
        %widefield post-processing parameters
        w_deconvolution = 'lucric'; %type of deconvolution (if none, then will filter by below parameters)
        w_filter_freq = [0.1 4]; %frequency range to filter widefield data
        w_filter_type = 'lowpass'; 
        w_normalization_method = 'full'; %pixelwise, full, or bounded
        w_norm_val = [0, 97.5]; %either the precentile or the value (if bounded) to normalize to
        w_chunk_dur = 300 %duration of training/testing chunks for fitting seqNMF in seconds
        
        
        %CNMF Defaults
        K = 25;              %10                                  %Number of factors
        L = 15;              %100                                 %Length (timebins) of each factor exemplar
        maxiter =300;        %100                                 %Maximum # iterations to run for seqNMf
        maxiter_fitlambda = 100; %50                               %max number of iterations for fitting when finding lambda (don't need to run to completion since relative) 
        lambda = 0;          %0.0005                              %Regularization parameter
        lambdaL1H = 0;       %0                                   %L1 sparsity parameter; Increase to make H's more sparse
        lambdaOrthoH = 1;    %0                                   %||HSH^T||_1,i~=j; Encourages events-based factorizations
        tolerance = 0;       %0                                   %Stop fitting if error reaches said value               
        shift = 1;    
        showPlot = 0; 
        lambda_range = sort(logspace(0,-6,15), 'ascend'); %in logspace range(1), range(2), number of lambdas
        
        %deconvolution defaults
        d_gamma = 0.975; 
        d_smooth = 1;
        d_kernel = 200;

    end

    methods
 
    end

end











