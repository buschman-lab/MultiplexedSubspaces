function [net,stats,netc,xic,aic,tr] = train_narx_nn(dff,fr,params)
%Camden MacDowell - timeless, adapted from auto matlab output'
%   autogressive network with exogenous input. 
%   dff - vector of calcium data. Reccomend to flip in time so predicting
%   firing rate from FUTURE imaging data. 
%   fr - vector of firing rates time series.
% returns the open loop network (net) and performance stats, closed loop (netc) and closed loop
% states (xic,aix)
if nargin <3
    params.n_hiddenlayer = 20; %neurons in hidden layer
    params.trainFcn = 'trainlm'; % 'trainlm' is usually fastest.
    params.inputdelay = 15; %number of timepoints for prediction
    params.feedbackdelay = 3; %number of timepoints for prediction
    params.verbose = 0; 
    params.hiddenfnc = 'tansig';
    params.outputfnc = 'purelin';
end

%List of transfer functions
% %     compet - Competitive transfer function.
% %     elliotsig - Elliot sigmoid transfer function.
% %     hardlim - Positive hard limit transfer function.
% %     hardlims - Symmetric hard limit transfer function.
% %     logsig - Logarithmic sigmoid transfer function.
% %     netinv - Inverse transfer function.
% %     poslin - Positive linear transfer function.
% %     purelin - Linear transfer function.
% %     radbas - Radial basis transfer function.
% %     radbasn - Radial basis normalized transfer function.
% %     satlin - Positive saturating linear transfer function.
% %     satlins - Symmetric saturating linear transfer function.
% %     softmax - Soft max transfer function.
% %     tansig - Symmetric sigmoid transfer function.
% %     tribas - Triangular basis transfer function.

X = tonndata(dff,true,false);
T = tonndata(fr,true,false);

% Create a Nonlinear Autoregressive Network with External Input
inputDelays = 1:params.inputdelay;
feedbackDelays = 1:params.feedbackdelay;
net = narxnet(inputDelays,feedbackDelays,params.n_hiddenlayer,'open',params.trainFcn);
net.layers{1}.transferFcn = params.hiddenfnc; %hidden layer
net.layers{2}.transferFcn = params.outputfnc; %hidden layer

% Choose Input and Feedback Pre/Post-Processing Functions
% Settings for feedback input are automatically applied to feedback output
% For a list of all processing functions type: help nnprocess
% Customize input parameters at: net.inputs{i}.processParam
% Customize output parameters at: net.outputs{i}.processParam
for i = 1:net.numInputs
    net.inputs{i}.processFcns = {'removeconstantrows'};
end
% net.inputs{2}.processFcns = {'removeconstantrows','mapminmax'};

% Prepare the Data for Training and Simulation
% The function PREPARETS prepares timeseries data for a particular network,
% shifting time by the minimum amount to fill input states and layer
% states. Using PREPARETS allows you to keep your original time series data
% unchanged, while easily customizing it for networks with differing
% numbers of delays, with open loop or closed loop feedback modes.
[x,xi,ai,t] = preparets(net,X,{},T);

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivision
net.divideFcn = 'divideblock';  % Divide data maintaining temporal relationships
net.divideMode = 'time';  % Divide up every sample
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0.2; %this is to prevent overfitting
net.divideParam.testRatio = 0.1; %doing this outside of this function

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
stats.train_indx = tr.trainMask;

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

%convert to closed loop
[netc,xic,aic] = closeloop(net,xf,af);
% %retrain
% x = preparets(netc,X,{},T);
% [netc, ~, ~, ~, xic, aic] = train(netc,x,t,xic,aic);
% %get performance
% [y,~,~] = netc(x,xic,aic);
% stats.performance = perform(netc,t,y);

end %function


