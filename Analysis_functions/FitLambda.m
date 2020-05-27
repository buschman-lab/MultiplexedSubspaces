function [lambda, fh] = FitLambda(X,gp,genfigs)
%Camden MacDowell - timeless. Follows Mackevivious et al., 2018, eLife
if nargin <3 
    genfigs =1;
end

fprintf('\n Fitting Lambda')
cost = NaN(1,numel(gp.lambda_range));
reg = NaN(1,numel(gp.lambda_range));
for i = 1:numel(gp.lambda_range)
    tic;
    fprintf('\n\t ........ Fitting Lambda %d of %d ...........', i,numel(gp.lambda_range));
    [w, h] = seqNMF(X, ...    
        'K', gp.K, 'L',gp.L, 'lambda',gp.lambda_range(i),...        
        'showPlot', 0, 'maxiter',gp.maxiter,'tolerance',gp.tolerance,...
        'SortFactors',0,'lambdaL1H',gp.lambdaL1H,...
        'lambdaOrthoH',gp.lambdaOrthoH,'useWupdate',1,'Shift',gp.shift);
    
    Xhat = tensor_convolve(w,h);
    [cost(i),reg(i),~] = helper.get_seqNMF_cost(X,w,h);
end

%normalize
cost_norm = (cost-min(cost))/(max(cost)-min(cost));
reg_norm = (reg-min(reg))/(max(reg)-min(reg));

%interpolate and find intersection (where crosses zero). could also do this by fitting a polynomial
xq = linspace(gp.lambda_range(1),gp.lambda_range(end),1000); 
cost_norm_int = interp1(gp.lambda_range,cost_norm,xq,'linear');
reg_norm_int = interp1(gp.lambda_range,reg_norm,xq,'linear');

%approximate intersection
idx = find(reg_norm_int-cost_norm_int<0,1,'first');
lambda = xq(idx);

%optionally make figures;
if genfigs
   figure; hold on; 
   plot(gp.lambda_range,cost_norm,'ro'); 
   plot(gp.lambda_range,reg_norm,'bo'); 
   plot(gp.lambda_range,pev_norm,'go'); 
   set(gca,'xscale','log');
   plot(xq,cost_norm_int,'r');plot(xq,reg_norm_int,'b'); plot(lambda,reg_norm_int(idx),'kx')
   fh = gcf;
else
   fh = [];
end

end
