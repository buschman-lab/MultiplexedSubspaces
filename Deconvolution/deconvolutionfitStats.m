function stats = deconvolutionfitStats(ypred,y)
%Camden MacDowell - timeless
%Computes basic stats about the quality of fit between the true and
%predcited

[n,z] = size(y);

%loop through all comparisons in ypred
stats = struct();
for cur_comp = 1:numel(ypred)    
    %correlation between predicted and true
    stats.rho(:,cur_comp) = diag(corr(ypred{cur_comp},y));

    %maximum crosscorrelation between pred and true (if accurately timed,
    %should match max_rho and lag=0;
    x_rho = NaN(z,1);
    x_lag = NaN(z,1);
    for i = 1:z
        [a,b] = xcorr(ypred{cur_comp}(:,i)-nanmean(ypred{cur_comp}(:,i)),y(:,i)-nanmean(y(:,i)),60,'normalized');
        [~,idx] = max(a);
        x_rho(i) = a(idx);
        x_lag(i) = b(idx);
    end
    stats.x_rho(:,cur_comp) = x_rho;
    stats.x_lag(:,cur_comp) = x_lag;

    %rmse between true and predicted
    err = NaN(z,1);
    for i = 1:z
        xhat = ypred{cur_comp}(:,i);
        x = y(:,i);
        err(i) = sqrt((1/n)*nansum(x-xhat)^2);
    end
    stats.err(:,cur_comp) = err;
    
end% ypred loop


end %function end