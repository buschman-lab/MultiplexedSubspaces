classdef behavioral_params  
    properties               
        %preprocessing parameters
        zscore = 1; 
        include_facecam = 0;
        derivative = 0;
        
        %ROI names. first cell is face movies, second is body movies. 
        roi_names = { {'timing','noise','eye'}, {'timing'} };  
        
        %Make behavioral dffs
        method = 'movingavg';
        detrend = 1; 
        fps = 60; 
        window  = 30; %in secconds
        
        %dlc parts list
        dlc_parts_list = {'nosetip','frontrightpawcenter','frontleftpawcenter',...
            'backrightpawcenter','backleftpawcenter','tail_1'};
        dlc_reference_part = {'tailroot'};
        dlc_epsilon = 50; 
                 
        %embedding parameters
        umap_mind_dist = 0.9;
        umap_n_neighbors = 50; 
        
        %Phenograph paramters
        pheno_knn = 50; 
        
        %motif trigger analysis parameters
        trig_dur = 10*13; %in frames post 
        min_std = 1;
        max_std = 2;
        pheno_level = 4;

        
        
    end

    
    
    
    methods
        
        
    end

end











