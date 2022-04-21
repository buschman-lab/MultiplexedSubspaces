function pev = PercentExplainedVariance(Ytest, Yhat, avgflag, lossFlag)
% Camden - timeless
%if avgflag =0 (def: computes across the entire population. If =1 then
%computes per neuron and averages

[n,K] = size(Ytest);

numModels = size(Yhat, 2)/K;

res = reshape( bsxfun(@minus, repmat(Ytest, [1 numModels]), Yhat), n,K,numModels);

pev = NaN(1,size(res,3));
for i = 1:size(res,3)
    if avgflag        
        temp = res(:,:,i); 
        pev(i) = nanmean(arrayfun(@(n) 1-nanvar(temp(:,n))./nanvar(Ytest(:,n)), 1:size(temp,2),'UniformOutput',1));        
    else
        %compute across the population
        temp = res(:,:,i); 
        pev(i) = 1-nanvar(temp(:))./nanvar(Ytest(:));
    end
end

if lossFlag %return as loss
    pev = 1-pev*100; 
else
    pev = pev*100; 
end

end