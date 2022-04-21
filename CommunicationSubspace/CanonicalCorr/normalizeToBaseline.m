function x_norm = normalizeToBaseline(x,base_tp,type)
%Camden MacDowell
%normalize trial-triggered spiking to average of baseline period
%x is a neuron x timepoint x trial tensor
%base_tp is the baseline timepoints

% base_tp = [1:5]

switch type
    case 'zscore'
        x_norm = (x-nanmean(x(:,base_tp,:),2))./(nanstd(x(:,base_tp,:),[],2)+1);
    case 'mean'
        x_norm = x./(nanmean(x(:,base_tp,:),2)+1);
    case 'meansubtract'
        x_norm = x-nanmean(x(:,base_tp,:),2);
    case 'none'
        x_norm = x;

end


end

