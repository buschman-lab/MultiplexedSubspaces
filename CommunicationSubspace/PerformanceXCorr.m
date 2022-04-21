function [xrho,idx,rho] = PerformanceXCorr(x,y,k,verbose)
%Camden MacDowell - timeless
% is the lag to use (k=1 (def) is the best, k=2 would be the second best...
% allows testing for dovetailing
if nargin <3; k = 1; end
if nargin <4; verbose = 0; end

[n,m,z] = size(x);
idx = NaN(n,m); 
xrho = NaN(n,m);
rho = NaN(n,m);
rho_trace_all = NaN(n,m,z*2-1);
for cur_rec = 1:n %rec loop
    for cur_m = 1:m %motif loop
        xx = squeeze(x(cur_rec,cur_m,:));
        xx = xx-nanmean(xx);
        yy = squeeze(y(cur_rec,cur_m,:));
        yy = yy-nanmean(yy);
        [rho_trace,lags] = xcorr(xx,yy,'normalized');
        [rho_temp,best_lag] = maxk(rho_trace,k);
        rho_temp = rho_temp(k); best_lag = best_lag(k);
        
        xrho(cur_rec,cur_m) = rho_temp;
        rho(cur_rec,cur_m) = corr(xx,yy);
        idx(cur_rec,cur_m) = lags(best_lag);
        rho_trace_all(cur_rec,cur_m,:) = rho_trace;
    end    
end

if verbose    
    figure; hold on; 
    for cur_m = 1:m
        nexttile
        plot(lags,squeeze(nanmean(rho_trace_all(:,cur_m,:),1)));
        title(sprintf('motif %d',cur_m),'fontweight','normal');
    end
end

%to visualize the lag
% figure; hold on; 
% plot(xx); 
% yyaxis right; plot(yy,'r')
% ytemp = circshift(yy,-1);
% plot(ytemp,'c'); 
end %function end