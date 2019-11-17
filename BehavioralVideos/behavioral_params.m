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
        dlc_epsilon = 15;                 
        
        %Make behavioral dffs
        method = 'movingavg';
        detrend = 1; 
        fps = 60; 
        window  = 30; %in secconds
        
        %motif triggered analysis parameters
        classification_idx = (1:10*13); %(1:26*30)
        trig_dur = (-13*5:13*5); %(-13*30:13*30); 
        min_std = 1;
        baseline = 1:2*13;
        noise_dur = (-45*13:-35*13);
        %@(x) randi(60*55)*13+(1:13*10+1) %select random points through recording
        %(randi(30)*13-(60*13))+(1:13*6) %select random range between -30 and -60 seconds
        
    end

    
    
    
    methods
        
        
    end

end











