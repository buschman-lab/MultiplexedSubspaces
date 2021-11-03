function [net,stats,tr] = train_feedforward_nn(dff,fr,params)
%Camden MacDowell - timeless 
if nargin <3
    params.n_hiddenlayer = 10; %neurons in hidden layer
    params.trainFcn = 'trainlm'; % 'trainlm' works best
    params.win = 60; %size of timepoints to use. 60 is best
    params.verbose = 0; 
    params.hiddenfnc = 'tansig'; %tansig, softmax, and radbasn, are all good
    params.outputfnc = 'purelin'; %pure linear performs best, poslin is comparable but forces postive output. Tansig is also okay
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

x = createRollingWindow(dff, params.win)'; %t-n:t-1
t =  fr(ceil(params.win/2):end-floor(params.win/2)); % get the middle timepoint in window  

% Create a Nonlinear Autoregressive Network with External Input
net = feedforwardnet(params.n_hiddenlayer,params.trainFcn);
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

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivision
net.divideFcn = 'divideblock';  % Divide data maintaining temporal relationships
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 0.8; %0.7
net.divideParam.valRatio = 0.2; %0.2 %this is to prevent overfitting
net.divideParam.testRatio = 0; %0.1 %doing this outside of this function

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate', 'ploterrhist', ...
    'plotregression', 'plotresponse', 'ploterrcorr', 'plotinerrcorr'};

% Train the Network
net.trainParam.showWindow = 0;
[net,tr] = train(net,x,t);

% Test the Network and get the open loop final states
y = net(x);
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

end %function


