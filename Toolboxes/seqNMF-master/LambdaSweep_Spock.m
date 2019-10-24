% function TestPrimitives
% Runs on spock
%Pipeline to take trained factors, test their signficance
%make plot for individual animal recs

function LambdaSweep_Spock(lambdablock)
%Load train discovered primitives
bucket = '/jukebox/buschman/';
base = 'Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/AnalyzedData/';
load([bucket base 'TrainRepitoires/13-Feb-2019_REPITOIRE_L3e4_H1.mat'],'H','W_all');
%Load the testing data
load([bucket base 'AllRecDFF/AllData_binned_SmallMask_3x_2minTesting.mat']);
%Load the significance information
load([bucket base 'TestRepitoires/2-15-2019-3x2min-chunks_Significance/' sprintf('SignificantFactors_chunk_%d.mat',block)]);
%Set Save Dir
savedir = [bucket base 'TestRepitoires/2-16-2019-3x2min-chunk_FitSignificantFactors/'];
%%
%Get the training chunk ID that corresponds to each factor in W_all (the number of factors for each chunk is contained in each cell in H)
Wid = cell(1,numel(H));
for cur_train = 1:numel(H)
    Wid{cur_train} = ones(1,size(H{cur_train},1))*cur_train;
end
Wid = [Wid{:}]';

%Recondition the test data
temp = conditionDffMat(data{block}',nanpxs{block});
data_test = reshape(temp,[size(temp,1)*size(temp,2),size(temp,3)]);

block = 10;

%get training and test data
data_train = W_all(:,Wid == block,:);

%Just get the significant factors
data_train = data_train(:,is_sig,:);

%Nan out pixels that have zero variance in test data
data_test(nanvar(data_test, [], 2) <= eps) = NaN;

%Remove all NaN pixel from both testing and training data
bad_pxl = isnan(data_test(:,1));
data_test(bad_pxl,:) = [];
data_train(bad_pxl,:,:) = [];

%Lambda Sweep 
[cost, reg, lambda] = ParamSweepFitted(data_train,data_test,lambdablock);


%Save off parameters

filename = [savedir sprintf('lambda_block_%d.mat',block)];
save(filename,'cost','reg','lambda','-v7.3');
end %function end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
