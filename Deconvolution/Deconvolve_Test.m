function [stats,stPred] = Deconvolve_Test(trace_probe,spikes_probe,method,trained_opts)
%Camden MacDowell - timeless
%applies deconvolution of the desired method on testing data
%trace_probe and spikes_probe are vectors of of dff and spikes of one probe 
timepoints = zeros(1,size(spikes_probe,1)); %track used timepoints (varies by method)
switch method
    case 'narx' %autoregressive neural network 
        params = trained_opts.narxparams;
        net =trained_opts.closedloopnetwork;
        X = tonndata(trace_probe',true,false);
        [x,xic,~,~] = preparets(net,X,{});
        aic = cat(1,repmat({zeros(params.n_hiddenlayer,1)},1,params.feedbackdelay),repmat({zeros(1,1)},1,params.feedbackdelay)); %init feedback condition
        stPred = net(x,xic,aic);        
        stPred = ([stPred{:}]');        
        timepoints(params.inputdelay+1:end)=1; %adjust for middle time points    
        
    case 'feedforward' %feedforward neural network   
        params = trained_opts.feedforwardparams;
        net =trained_opts.shallowfeedforward;
        x = createRollingWindow(trace_probe', params.win)'; %t-n:t-1        
        stPred = net(x)';        
        timepoints(ceil(params.win/2):end-floor(params.win/2))=1; % get the middle timepoint in window  
    case 'lr_glm' %lucy richarson deconvolution with glm kernel
        kernel = trained_opts.glmkernel; 
        kernel = kernel-min(kernel); %requires positive
        stPred = deconvlucy(trace_probe-min(trace_probe),kernel); 
        timepoints(1:end)=1;
    case 'lr_gcamp' %lucy richarson deconvolution with gcamp kernel
        win = trained_opts.LRwin;
        gamma = trained_opts.LRgamma;
        stPred = lucric(trace_probe-min(trace_probe),gamma,1,win); 
        timepoints(1:end)=1;         
    case 'glm'
        kernel = flipud(trained_opts.glmkernel);
        stPred = convn(padarray(trace_probe',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid');  
        timepoints(1:end)=1;
    case 'none'
        stPred = trace_probe;
        timepoints(1:end)=1;        
    otherwise
        error('unknown deconvolv method')  
end %method switch

stats = deconvolutionfitStats(stPred,spikes_probe(timepoints==1));


end %function 
