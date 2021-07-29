function [net,stats,netc,xic,aic] = train_narx_nn(dff,fr,params)
%Camden MacDowell - timeless, adapted from auto matlab output'
%   autogressive network with exogenous input. 
%   dff - vector of calcium data. Reccomend to flip in time so predicting
%   firing rate from FUTURE imaging data. 
%   fr - vector of firing rates time series.
% returns the open loop network (net) and performance stats, closed loop (netc) and closed loop
% states (xic,aix)
if nargin <3
    params.n_hiddenlayer = 10; %neurons in hidden layer
    params.trainFcn = 'trainlm'; % 'trainlm' is usually fastest.
    % 'trainbr' takes longer but may be better for challenging problems.
    % 'trainscg' uses less memory. Suitable in low memory situations.
    params.delay = 20; %number of timepoints for prediction
    params.verbose = 0; 
end

X = tonndata(dff,true,false);
T = tonndata(fr,true,false);

% Create a Nonlinear Autoregressive Network with External Input
inputDelays = 1:params.delay;
feedbackDelays = 1:params.delay;
net = narxnet(inputDelays,feedbackDelays,params.n_hiddenlayer,'open',params.trainFcn);

% Choose Input and Feedback Pre/Post-Processing Functions
% Settings for feedback input are automatically applied to feedback output
% For a list of all processing functions type: help nnprocess
% Customize input parameters at: net.inputs{i}.processParam
% Customize output parameters at: net.outputs{i}.processParam
net.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
net.inputs{2}.processFcns = {'removeconstantrows','mapminmax'};

% Prepare the Data for Training and Simulation
% The function PREPARETS prepares timeseries data for a particular network,
% shifting time by the minimum amount to fill input states and layer
% states. Using PREPARETS allows you to keep your original time series data
% unchanged, while easily customizing it for networks with differing
% numbers of delays, with open loop or closed loop feedback modes.
[x,xi,ai,t] = preparets(net,X,{},T);

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivision
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'time';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 30/100; %this is to prevent overfitting
net.divideParam.testRatio = 0; %doing this outside of this function

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate', 'ploterrhist', ...
    'plotregression', 'plotresponse', 'ploterrcorr', 'plotinerrcorr'};

% Train the Network
net.trainParam.showWindow = 0;
[net,tr] = train(net,x,t,xi,ai);

% Test the Network and get the open loop final states
[y,xf,af] = net(x,xi,ai);
e = gsubtract(t,y);
stats.performance = perform(net,t,y);

% Recalculate Training, Validation and Test Performance
trainTargets = gmultiply(t,tr.trainMask);
valTargets = gmultiply(t,tr.valMask);
testTargets = gmultiply(t,tr.testMask);
stats.trainPerformance = perform(net,trainTargets,y);
stats.valPerformance = perform(net,valTargets,y);
stats.testPerformance = perform(net,testTargets,y);

% Plots
% Uncomment these lines to enable various plots.
if params.verbose
    view(net)
    figure, plotperform(tr)
    figure, plottrainstate(tr)
    figure, ploterrhist(e)
    figure, plotregression(t,y)
    figure, plotresponse(t,y)
    figure, ploterrcorr(e)
    figure, plotinerrcorr(x,e)
end

%simulate the closed loop network which can be run without spike input
[netc,xic,aic] = closeloop(net,xf,af);

end %function


