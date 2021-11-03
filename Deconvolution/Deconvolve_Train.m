function trained_opts = Deconvolve_Train(dff,st,method,min_spikes,bindata)
%Camden MacDowell - timeless
%creates a model or kernel for given deconovlution method
%use in tandem with Deconvolve_Test
%dff and st are cell arrays of dff and spikes per recording

if nargin <4
    min_spikes =2000; 
end
n_rec = numel(dff);

%adjust to 2 second windows for everything
if bindata ==1
    win=30;
else
    win=60;
end


trained_opts = cell(1,n_rec); %cell array of structures for each rec
%loop through recording
for cur_rec = 1:n_rec
    trained_opts_cur_rec = struct();
    trace = dff{cur_rec};
    spikes = st{cur_rec};       
    n_probe = size(trace,2);
    for cur_probe = 1:n_probe
        fprintf('\n\t training models probe %d of %d rec %d of %d',cur_probe,n_probe,cur_rec, n_rec)        
        %break out variables for ease
        trace_probe = trace(:,cur_probe); 
        spikes_probe = spikes(:,cur_probe); 
        if sum(spikes_probe) > min_spikes %make sure at least some spikeing activity
            switch method
                case 'narx' %autoregressive neural network                
                    %train
                    params = [];
                    params.n_hiddenlayer = 20; %neurons in hidden layer
                    params.trainFcn = 'trainlm'; % 'trainlm' is usually fastest.
                    params.inputdelay = win; %number of timepoints for prediction
                    params.feedbackdelay = 3; %number of timepoints for prediction
                    params.verbose = 0; 
                    params.hiddenfnc = 'tansig';
                    params.outputfnc = 'purelin';
                    [~,~,net,~,~,~] = train_narx_nn(trace_probe',spikes_probe',params);
                    %save off needed variables
                    trained_opts_cur_rec(cur_probe).closedloopnetwork = net;
                    trained_opts_cur_rec(cur_probe).narxparams = params;
                    %to test: stPred{i}(:,j) = cell2mat(netc(dff_test,xic,aic));
                case 'glm' %linear deconvolution using glm kernel ... this is least squares... can test by comparing kernel of fitlm with glmfit
                    winglm = win+1; %size of the kernel window in frames +1(def = 1 sec=30)
                    predictors = createRollingWindow(trace_probe, winglm); %t-n:t-1
                    response =  spikes_probe(ceil(winglm/2):end-floor(winglm/2)); % get the middle timepoint in window  
                    kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
                    trained_opts_cur_rec(cur_probe).glmkernel = kernel(2:end); %remove intercept  
                    trained_opts_cur_rec(cur_probe).glmkernelintercept = kernel(1);
                    %to test: stPred = {convn(padarray(dff',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid')'};   
                case 'feedforward' %linear deconvolution using glm kernel
                    %shallow feedforward
                    fprintf('\n\t feedforward')
                    params = [];
                    params.n_hiddenlayer = 20; %neurons in hidden layer
                    params.trainFcn = 'trainlm'; % 'trainlm' works best
                    params.win = win; %size of timepoints to use. 60 is best
                    params.verbose = 0; 
                    params.hiddenfnc = 'tansig'; %tansig, softmax, and radbasn, are all good
                    params.outputfnc = 'purelin'; %pure linear performs best, poslin is comparable but forces postive outpu               
                    net = train_feedforward_nn(trace_probe',spikes_probe',params); %train 
                    trained_opts_cur_rec(cur_probe).shallowfeedforward = net;
                    trained_opts_cur_rec(cur_probe).feedforwardparams = params; 
                case 'lr_gcamp' %adapt for temporal binning
                    if bindata ==1 %0.95 if 15fps, 0.96 if 30
                        trained_opts_cur_rec(cur_probe).LRgamma = 0.80; %predicted is 0.95
                    else
                        trained_opts_cur_rec(cur_probe).LRgamma = 0.89; %predicted is 0.96
                    end
                    trained_opts_cur_rec(cur_probe).LRwin = win;
                case 'all'
%                     %narx
%                     fprintf('\n\t narx')
%                     params = [];
%                     params.n_hiddenlayer = 20; %neurons in hidden layer
%                     params.trainFcn = 'trainlm'; % 'trainlm' is usually fastest.
%                     params.inputdelay = win; %number of timepoints for prediction
%                     params.feedbackdelay = 3; %number of timepoints for prediction
%                     params.verbose = 0; 
%                     params.hiddenfnc = 'tansig';
%                     params.outputfnc = 'purelin';
%                     [~,~,net,~,~,~] = train_narx_nn(trace_probe',spikes_probe',params);
%                     %save off needed variables
%                     trained_opts_cur_rec(cur_probe).closedloopnetwork = net;
%                     trained_opts_cur_rec(cur_probe).narxparams = params;

                    %glm
                    fprintf('\n\t glm')
                    winglm = win+1; %size of the kernel window in frames (def = 1 sec=30)
                    predictors = createRollingWindow(trace_probe, winglm); %t-n:t-1
                    response =  spikes_probe(ceil(winglm/2):end-floor(winglm/2)); % get the middle timepoint in window  
                    kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
                    trained_opts_cur_rec(cur_probe).glmkernel = kernel(2:end); %remove intercept  
                    trained_opts_cur_rec(cur_probe).glmkernelintercept = kernel(1);

                    %shallow feedforward
                    fprintf('\n\t feedforward')
                    params = [];
                    params.n_hiddenlayer = 20; %neurons in hidden layer
                    params.trainFcn = 'trainlm'; % 'trainlm' works best
                    params.win = win; %size of timepoints to use. 60 is best
                    params.verbose = 0; 
                    params.hiddenfnc = 'tansig'; %tansig, softmax, and radbasn, are all good
                    params.outputfnc = 'purelin'; %pure linear performs best, poslin is comparable but forces postive outpu               
                    net = train_feedforward_nn(trace_probe',spikes_probe',params); %train 
                    trained_opts_cur_rec(cur_probe).shallowfeedforward = net;
                    trained_opts_cur_rec(cur_probe).feedforwardparams = params;

                    %lucrid
                    if bindata ==1 %0.95 if 15fps, 0.96 if 30
                        trained_opts_cur_rec(cur_probe).LRgamma = 0.80; %predicted is 0.95
                    else
                        trained_opts_cur_rec(cur_probe).LRgamma = 0.89; %predicted is 0.96
                    end
                    trained_opts_cur_rec(cur_probe).LRwin = win;
                otherwise 
                    error('unknown method');
            end %method switch
            trained_opts_cur_rec(cur_probe).insufactivity = 0;
        else
            trained_opts_cur_rec(cur_probe).insufactivity = 1;
        end %number of spikes            
    end %probe loop
    trained_opts{cur_rec} = trained_opts_cur_rec;
end %rec loop


end %function 
