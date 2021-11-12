function stats = deconvolutionfitStats(ypred,y)
%Camden MacDowell - timeless
%Computes basic stats about the quality of fit between the true and
%predcited
%y and ypred is a vector 

n = size(y,1);
  
%correlation between predicted and true (diag = match probes)
stats.rho = diag(corr(ypred-nanmean(ypred),y-nanmean(y)));

%maximum crosscorrelation between pred and true 
[a,b] = xcorr(ypred-nanmean(ypred),y-nanmean(y),30,'normalized');
[~,idx] = max(a);
stats.x_rho = a(idx);
stats.x_lag = b(idx);

%get a long window for plotting
[a,~] = xcorr(ypred-nanmean(ypred),y-nanmean(y),600,'normalized');
stats.xcorrvect = a;
[a,~] = xcorr(ypred-nanmean(ypred),600,'normalized');
stats.autocorpred = a;
[a,~] = xcorr(y-nanmean(y),600,'normalized');
stats.autocortrue = a;

%compute the skewness as in https://en.wikipedia.org/wiki/Skewness 
%so the third moment of x/(second moment .^1.5). 
%will be negative if neg skewed, pos if pos skewed and 0 if no skew
x = ypred-nanmean(ypred); %center
stats.skew_pred = nanmean(x.^3) ./ (nanmean(x.^2).^1.5);

x = y-nanmean(y); %center
stats.skew_true = nanmean(x.^3) ./ (nanmean(x.^2).^1.5);

%difference in skew
stats.s = stats.skew_true-stats.skew_pred;

%rmse between true and predicted
stats.err = sqrt((1/n)*nansum((y-ypred).^2));

end %function end