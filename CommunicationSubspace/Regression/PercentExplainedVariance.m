function pev = PercentExplainedVariance(Ytest, Yhat)
% 

[~, K] = size(Ytest);

numModels = size(Yhat, 2)/K;

pev = NaN(1,K);
for i = 1:K
    resid = Ytest-Yhat(:,i);
    pev(i) = 1 - nanvar(resid(:))./nanvar(Ytest(:));
end

end