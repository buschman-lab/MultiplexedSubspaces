function spock_compareNNparams(type,block)
fprintf('working on block %d',block);

%addpath
addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));

%load data
% temp = load('/jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/analysisplayground/temp_data/decondata.mat');
temp = load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/analysisplayground/temp_data/decondata.mat');
trace_probe_train=temp.trace_probe_train;
spikes_probe_train = temp.spikes_probe_train;


switch type
    case 1
        fprintf('NARX')
        %CreateDifferentData
        n_hiddenlayer = [{20},{5,10},{5,20},{4,5}];
        inputdelay = [15,30,60,120]; 
        feedbackdelay = [2,4,6,20]; 
        trainFcn = {'trainlm'}; % 'trainbr'
        hiddenfunc_list = {'elliotsig','softmax','tansig'};
        outputfunc_list = {'softmax','purelin'};
        % trainFcn = {'trainlm','trainbr'}; % 'trainbr'
        % hiddenfunc_list = {'elliotsig','hardlim','poslin','radbas','radbasn','satlin','softmax','tansig','tribas','purelin'};
        % outputfunc_list = {'elliotsig','hardlim','hardlims','logsig','netinv','poslin','radbas','radbasn','satlin','softmax','tansig','tribas','purelin'};
        totest = combvec(1:numel(n_hiddenlayer),1:numel(trainFcn),inputdelay,feedbackdelay,1:numel(hiddenfunc_list),1:numel(outputfunc_list));

        fprintf('working on block %d of %d',block,size(totest,2));

        cur_test = totest(:,block);

        %load the example data
        params.n_hiddenlayer = n_hiddenlayer{cur_test(1)}; %neurons in hidden layer
        params.trainFcn = trainFcn{cur_test(2)};
        params.inputdelay = cur_test(3); %number of timepoints for prediction
        params.feedbackdelay = cur_test(4); %number of timepoints for predictionparams.verbose = 0; 
        params.hiddenfnc = hiddenfunc_list{cur_test(5)};
        params.outputfnc = outputfunc_list{cur_test(6)};
        params.verbose = 0; 

        % fit and test model

        rng('default');
        X = tonndata(trace_probe_train',true,false);
        fprintf('\n\t working');
        [~,~,netc,~,~,tr] = train_narx_nn(trace_probe_train',spikes_probe_train',params);
        [x,xic,~,~] = preparets(netc,X,{});
        aic = cat(1,repmat({zeros(params.n_hiddenlayer,1)},1,params.feedbackdelay),repmat({zeros(1,1)},1,params.feedbackdelay)); %init feedback condition
        y_train = netc(x,xic,aic);        
        y_train = ([y_train{:}]');
        temp = spikes_probe_train(params.inputdelay+1:end);
        t_true = temp(1:numel(y_train));
        rho = corr(y_train(tr.trainInd)',t_true(tr.trainInd)');
        e = mse(y_train(tr.trainInd)',t_true(tr.trainInd)');     
        rho_gen = corr(y_train([tr.valInd,tr.testInd])',t_true([tr.valInd,tr.testInd])');
        e_gen = mse(y_train([tr.valInd,tr.testInd])',t_true([tr.valInd,tr.testInd])');    

    case 2 %feedfoward
        %CreateDifferentData
        fprintf('FeedForward')
        n_hiddenlayer = [{20},{10},{40},{80},{100}];
        win = [60];
        trainFcn = {'trainlm'}; % 'trainbr'
        hiddenfunc_list = {'radbasn','softmax','tansig'};
        outputfunc_list = {'poslin','purelin','tansig'};
        % trainFcn = {'trainlm','trainbr'}; % 'trainbr'
%         hiddenfunc_list = {'poslin','radbas','radbasn','satlin','softmax','tansig','tribas','purelin'};
%         outputfunc_list = {'hardlims','logsig','netinv','poslin','radbas','radbasn','satlin','softmax','tansig','tribas','purelin'};

        totest = combvec(1:numel(n_hiddenlayer),1:numel(trainFcn),win,1:numel(hiddenfunc_list),1:numel(outputfunc_list));

        fprintf('working on block %d of %d',block,size(totest,2));

        cur_test = totest(:,block);

        %load the example data
        params.n_hiddenlayer = n_hiddenlayer{cur_test(1)}; %neurons in hidden layer
        params.trainFcn = trainFcn{cur_test(2)};
        params.win = cur_test(3);
        params.hiddenfnc = hiddenfunc_list{cur_test(4)};
        params.outputfnc = outputfunc_list{cur_test(5)};
        params.verbose = 0; 
        
        %fit and test model
        rng('default');
        [netc,~,tr] = train_feedforward_nn(trace_probe_train',spikes_probe_train',params); %train 
        fprintf('\n\tdone training')
        x = createRollingWindow(trace_probe_train', params.win)'; %t-n:t-1
        t_true =  spikes_probe_train(ceil(params.win/2):end-floor(params.win/2))'; % get the middle timepoint in window  
        y_train = netc(x);
        rho = corr(y_train(tr.trainInd)',t_true(tr.trainInd)');
        e = mse(y_train(tr.trainInd)',t_true(tr.trainInd)');     
        rho_gen = corr(y_train([tr.valInd,tr.testInd])',t_true([tr.valInd,tr.testInd])');
        e_gen = mse(y_train([tr.valInd,tr.testInd])',t_true([tr.valInd,tr.testInd])');             
        

    otherwise 
        error('unknown type');
end

fprintf('\n\t saving')        
save(sprintf('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/analysisplayground/temp_data/block%d.mat',block),'rho','e','y_train','t_true','params','block','rho_gen','e_gen')
fprintf('\n\t done')


%%NONSPOCKDATA
% dff_list = GrabFiles('dff_combined.mat',1,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging'});
% spike_opts_list = GrabFiles('ap_opts.mat',1,{'H:\'});
% Load example data
% load the data
% params.mua = 1; %1= use both 'good' and 'mua' units. 0 = just 'good'
% params.depth = [0 250]; %depth from surface of probe
% params.radius = 2; %pixel radius around probe tip    
% params.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]
% params.splittype = 'half'; %type of training/test split to use (current just first and second half of recording
% [dff,st,~] = CompileData_deconvolution(dff_list(4),spike_opts_list(4),params);
% %split data 
% cur_probe = 1; 
% n = size(dff{1},1)*3;
% trace_probe_train = dff{1}(1:floor(n/4),cur_probe); 
% spikes_probe_train = st{1}(1:floor(n/4),cur_probe); 
% save('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\temp_data\decondata.mat','trace_probe_train','spikes_probe_train');
end %function end