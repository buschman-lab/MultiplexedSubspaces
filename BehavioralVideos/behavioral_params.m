classdef behavioral_params  
    properties               
        %preprocessing parameters
        zscore = 1; 
        include_facecam = 1;
        derivative = 1;
        
        %ROI names. first cell is face movies, second is body movies. 
        roi_names = { {'timing','nose','whiskerpad'}, {'timing'} };         
        
        %dlc parts list
        dlc_parts_list = {'nosetip','frontrightpawcenter','frontleftpawcenter',...
            'backrightpawcenter','backleftpawcenter','tailroot'};
        dlc_reference_part = [];
        dlc_epsilon = 50;                 
        
        %motif trigger analysis parameters
        trig_dur = (-13*3:13*3); 
        min_std = 1;
        baseline = 1:2*13;
        
        
    end

    
    
    
    methods
        
        
    end

end











