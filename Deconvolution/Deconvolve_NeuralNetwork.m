function [stPred,stTrue,stats] = Deconvolve_NeuralNetwork(dff,st,train_idx)
%Camden MacDowell - timeless
%set train index to 1 to train on everything

[n,z] = size(dff);

if nargin <3
    train_idx = false(n,1);
    train_idx(1:floor(n*0.5),:) = true;
end

%flip timecourse since we are predicting 'future'
dff = flipud(dff);
st = flipud(st);

if isempty(train_idx) %no comparisons
    dff = cat(1,dff(:));
    st = cat(1,st(:));    
    [~,~,netc,xic,aic] = train_narx_nn(dff',st');
    stPred = cell2mat(netc(dff',xic,aic));                        
    stTrue = st';
    
    %reflip and pad final 20 timepoints (size of predictive kernel)
    k = z-size(stPred,1);
    stTrue = [flipud(stTrue),NaN(k,1)];
    stPred = {[flipud(stPred),NaN(k,1)]};  
    
else %run comparisons
    stPred = cell(1,z);
    for i = 1:z %loop through probes/recs
        fprintf('\n\t running narx fits for probe/rec %d of %d', i,z);
        %train network        
        [~,stats,netc,xic,aic] = train_narx_nn(dff(train_idx,i)',st(train_idx,i)');
        
        stPred{i} = NaN(sum(train_idx==0),z);
        %test on all (including self)
        for j = 1:z %interal loop through probes/recs
            dff_test = num2cell(dff(~train_idx,j)');
            stPred{i}(:,j) = cell2mat(netc(dff_test,xic,aic));                        
        end %j loop
        stPred{i} = flipud(stPred{i});
    end %i loop
    stTrue = flipud(st(~train_idx,:));
    
end %comparison if/else



end %function 

   
%    %NNT | reverse so predicting st from imaging
%    win = floor(numel(st)/2);
%    [net,stats,netc,xic,aic] = train_narx_nn(flipud(dff_probe(1:win,cur_probe))',flipud(st(1:win))');
%    
%    %predict the second half of the recording using closed loop autoregression
%    testdata = num2cell(flipud(dff_probe(win+1:end,cur_probe))');
%    actdata = num2cell(flipud(st(win+1:end))');
%    yPred = netc(testdata,xic,aic);
%    
%    figure; plot(cell2mat(yPred),'color','k'); hold on; plot(cell2mat(actdata),'color','r')
%    figure; hold on; plot(xcorr(cell2mat(yPred)',cell2mat(actdata)',1000,'normalized'))
%   
%    yPredtrain = netc(num2cell(flipud(dff_probe(1:win,cur_probe))'),xic,aic);
%    actdatatrain = num2cell(flipud(st(1:win))');
%    figure; plot(cell2mat(yPredtrain),'color','k'); hold on; plot(cell2mat(actdatatrain),'color','r')
%    figure; hold on; plot(xcorr(cell2mat(yPredtrain)',cell2mat(actdatatrain)',1000,'normalized'))   
%    
%    