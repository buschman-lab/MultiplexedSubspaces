function stats = CNMFStats(w,h,X,rmv_flag)
if nargin <4
    rmv_flag = 1; %remove nearly empty ws
end

%generate stats for a fit to the data
if rmv_flag %Remove empty w or w with barely anything
    indempty = sum(sum(w>0,1),3)==0; % W is literally empty
    Wflat = sum(w,3); 
    indempty = indempty | (max(Wflat,[],1).^2> .5*sum(Wflat.^2,1)); % or one pixel has >50% of the power
    w(:,indempty,:) = []; % Delete factors that meet the above critera
    h(indempty,:) = [];
end

stats.number_motifs = size(w,2);
Xhat = helper.reconstruct(w,h);
Residuals = X-Xhat;

%explained variance
stats.pev = CalculateExplainedVariance(X,Residuals);

%per frame
stats.pev_frame = CalculateExplainedVarianceFrameWise(X,Residuals);

%correlation
stats.rho = corr(X(:),Xhat(:));

%Correlation per frame
stats.rho_frame = zeros(1,size(Xhat,2));
for i = 1:size(Xhat,2)    
    stats.rho_frame(i) = corr(Xhat(:,i),X(:,i));
end    

%loadings of each w
K = size(h,1); 
stats.loadings = zeros(1,K);        
for i = 1:K
    temp = helper.reconstruct(w(:,i,:),h(i,:)); 
    stats.loadings(i) = CalculateExplainedVariance(X,X-temp);
end
stats.loadings = stats.loadings/sum(stats.loadings);

%final cost
stats.rmse = {sqrt(mean((X(:)-Xhat(:)).^2))};
    
    
    
end
    
    
    
    
    
    
    
    
    
    
