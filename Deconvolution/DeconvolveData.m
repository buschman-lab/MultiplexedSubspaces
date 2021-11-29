function ypred = DeconvolveData(y_all,method,trained_opts,addGLMint)
%Camden MacDowell - timeless
%accepts y_all as a vector or y_all as a t x n matrix
%some methods return a shorter vector, so also return the valid timepoints. 
if nargin <4 
    addGLMint = 0; 
end
ypred = NaN(size(y_all));
n = size(y_all,2);
for cur_trace = 1:n %trace loop
    timepoints = NaN(1,size(y_all,1)); %track used timepoints (varies by method)
    y = y_all(:,cur_trace);
    switch method %method switch
        case 'narx' %autoregressive neural network 
            params = trained_opts.narxparams;
            net =trained_opts.closedloopnetwork;
            X = tonndata(y',true,false);
            [x,xic,~,~] = preparets(net,X,{});
            aic = cat(1,repmat({zeros(params.n_hiddenlayer,1)},1,params.feedbackdelay),repmat({zeros(1,1)},1,params.feedbackdelay)); %init feedback condition
            ypred_temp = net(x,xic,aic);        
            ypred_temp = ([ypred_temp{:}]');        
            timepoints(params.inputdelay+1:end)=1; %adjust for middle time points            
        case 'feedforward' %feedforward neural network   
            params = trained_opts.feedforwardparams;
            net =trained_opts.shallowfeedforward;
            x = createRollingWindow(y', params.win)'; %t-n:t-1        
            ypred_temp = net(x)';        
            ypred_temp(ypred_temp<0)=0; %convert purelin to relu/poslin
            timepoints(ceil(params.win/2):end-floor(params.win/2))=1; % get the middle timepoint in window  
        case 'lr_glm' %lucy richarson deconvolution with glm kernel
            kernel = trained_opts.glmkernel; 
            kernel = kernel-min(kernel); %requires positive
            ypred_temp = deconvlucy(y-min(y),kernel); 
            timepoints(1:end)=1;
        case 'lr_gcamp' %lucy richarson deconvolution with gcamp kernel
            win = trained_opts.LRwin;
            gamma = trained_opts.LRgamma;
            ypred_temp = lucric(y-min(y),gamma,1,win); 
            timepoints(1:end)=1;           
        case 'glm'            
            fprintf('\naddGLMint %d',addGLMint);        
            params = trained_opts.feedforwardparams;
            x = createRollingWindow(y', params.win); %t-n:t-1  
            mdl = trained_opts.glmmodel;
            ypred_temp = predict(mdl,x);
            ypred_temp = ypred_temp+abs(min(ypred_temp)); %make all positive
            timepoints(ceil(params.win/2):end-floor(params.win/2))=1; % get the middle timepoint in window  
        case 'none'
            ypred_temp = y;
            timepoints(1:end)=1;        
        otherwise
            error('unknown deconvolv method')          
    end %method switch
    %pad with nan to match the length of original (for later chopping)
    timepoints(timepoints==1)=ypred_temp;
    ypred(:,cur_trace) = timepoints;
end %trace loop


end %function endfprintf('\naddGLMint %d',addGLMint);        

