classdef behavioral_params  
    properties               
        %preprocessing parameters
        zscore = 0; 
        include_facecam = 1;
        derivative = 1;
        
        %ROI names. first cell is face movies, second is body movies. 
        roi_names = { {'timing','nose','whiskerpad'}, {'timing'} };         
        
        %dlc parts list
        dlc_parts_list = {'nosetip','frontrightpawcenter','frontleftpawcenter',...
            'backrightpawcenter','backleftpawcenter','tailroot'};
        dlc_epsilon = [];     
        
        movmode_dur = 13; %movmode window (actually x2 since on both sided)
        
        %Make behavioral dffs
        method = 'movingavg';
        detrend = 1; 
        fps = 60; 
        window  = 30; %in secconds
        
        %motif triggered analysis parameters
        classification_idx = ((4*13)+1:(8*13)+1); %(1:26*30)
        plotting_idx =  ((1*13):(11*13));
        trig_dur = ((-13*1):13*1); %(-13*30:13*30); 
        min_std = 1;
        baseline = 1:2*13;
        noise_dur = (-49*13:-41*13);
        %@(x) randi(60*55)*13+(1:13*10+1) %select random points through recording
        %(randi(30)*13-(60*13))+(1:13*6) %select random range between -30 and -60 seconds
        
    end

    
    
    
    methods
        
        
    end

end











