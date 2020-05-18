classdef general_params  
    properties
        %Path Options
        local_bucket = 'Z:\';
        spock_bucket = '\jukebox\buschman\';
        repo_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\';        
        dynamic_script_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\Spock\DynamicScripts\';
        processing_intermediates = 'Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Scratch_Processing_Intermediates\'; %location of the intermediate files in the processing pipeline
        figure_save_directory = 'Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\ProcessingFigures\'; %directory to save off figures during processing
        
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
        denoise_powerfrac = 0.5;
        denoise_pxlfrac = 0.01; %fraction of pixels to use. 1% seems good. Can gut check with ImpactOfPCADenoise function in gutcheck directory
        w_deconvolution = 'lucric'; %type of deconvolution (if none, then will filter by below parameters)
        w_filter_freq = [0.1 5]; %frequency range to filter widefield data if using no deconvolution
        w_filter_type = 'lowpass'; 
        w_normalization_method = 'full'; %pixelwise, full, or bounded
        w_norm_val = 95; %either the precentile or the value (if bounded) to normalize to
        w_chunk_dur = 150 %duration of training/testing chunks for fitting seqNMF in seconds
        w_approx_chunk_num = ceil(51926/(150*15)/2); %(total duration/w_chunk_dur*fps)/2 (for test and train split) This is used in pipeline to parallelize motifs fittings spock jobs without knowing the exact chunk number used. Unused will just fail as spock jobs. 
        
        
        %CNMF Defaults
        K = 25;              %10                                  %Number of factors
        L = 15;              %100                                 %Length (timebins) of each factor exemplar
        maxiter =300;        %100                                 %Maximum # iterations to run for seqNMf        
        lambda = -1;          %0.0005                              %Regularization parameter. Set ot -1 to have autofit per run
        lambdaL1H = 0;       %0                                   %L1 sparsity parameter; Increase to make H's more sparse
        lambdaOrthoH = 1;    %0                                   %||HSH^T||_1,i~=j; Encourages events-based factorizations
        tolerance = 0;       %0                                   %Stop fitting if error reaches said value               
        shift = 0;    
        showPlot = 0; 
        lambda_range = sort(logspace(0,-5,15), 'ascend'); %in logspace range(1), range(2), number of lambdas
        
        %deconvolution defaults
        d_gamma = 0.975; 
        d_smooth = 2;
        d_kernel = 200;
        
        %clustering defaults
        m_smooth_kernel = []; %[X, Y, Z] spatio-temporal gaussian smoothing kernel for clustering motifs; leave empty for no smoothing. 
        m_maxshift = 7; %maximum temporal shift (in + and - directions) for clustering. def= floor(0.5*motif length);
        m_community_fraction = 0.50; %percent most interconnected motifs to used for creating basis motifs
        m_removepad = 0; %totally remove pad from temporally allignment (1 = results in motifs of same length as L, otherwise can be a little longer)
        pheno_louvain_restarts = 5; %number of random initializations of the louvain community detection during clustering. Takes the initialization with the best clustering (ovq_q) value. 
        pheno_k_range = [2,4,6,8,10,15,20,30,50]; %range of phenograph values to sweep
        pheno_num_resamples = 50; %number of random samples (with resampling) to take from the population when finding k
        
        %miscellaneous additions
        pixel_dim = [68,38];

    end

    methods
 
    end

end











