function [w, cost, loadings, power, numFactors, Xhat, Residuals, W, H] = ...
    CollectNMFOutputs(w, h, cost, loadings, power, nanpxs, block)

%Gather in case used GPU
cost = {gather(cost)};
loadings = {gather(loadings)};
power = {gather(power)};
w = gather(w);

%Remove empty w or w with barely anything
indempty = sum(sum(w>0,1),3)==0; % W is literally empty
Wflat = sum(w,3); 
indempty = indempty | (max(Wflat,[],1).^2> .5*sum(Wflat.^2,1)); % or one pixel has >50% of the power
w(:,indempty,:) = []; % Delete factors that meet the above critera
H = gather(h);
H(indempty,:) = [];

%Number of factors
numFactors = {size(w,2)};

%Get Residuals: 
Xhat = helper.reconstruct(w,H);
Residuals = X-Xhat; 

%Recondition W with the Nan pixels to allow spatial comparisons
W = zeros((size(w,1)+size(nanpxs{block},2)),size(w,2),size(w,3));
for cur_fact = 1:size(w,2)
       wtemp = squeeze(w(:,cur_fact,:));
       wtemp = conditionDffMat(wtemp',nanpxs{block}); 
       wtemp(isnan(wtemp))=0;
       W(:,cur_fact,:) = reshape(wtemp,[size(wtemp,1)^2,1,size(wtemp,3)]);
end

end %function end