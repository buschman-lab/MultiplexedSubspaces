function [metric,tau_hl,rmv_metric,rmv_tau] = DynamicMetric(w)

%Camden MacDowell 2019
%Dynamic metric computes spatail correlation between all timepoints in each
%motif W is a N x Frame x K tensor from seqNMF
%Caluclates how 'Dynamic' the input W is by taking the sum of the
%off-diagnol of the frame x frame correlation matrix

%Build a lag matrix 
lags = [(size(w,3):-1:2) (1:1:size(w,3))];
lagmat = zeros(1,size(w,3));
for i = 1:size(w,3)
%     lagmat(i,:) = lags(13:end);
    lagmat(i,:) = lags(size(w,3):end);
    lags = circshift(lags,1);    
end
    
%Calculate Spatial Correlation
metric = NaN(size(w,2),length(unique(lagmat)));
tau_hl = NaN(size(w,2),1);
for k = 1:size(w,2)
    data = squeeze(w(:,k,:));
    %Remove pixels that were not active in the motif (e.g. zero variance)
    bad_dim = nanvar(data,[],2)<=eps;    
    data(bad_dim,:) = []; %remove blank pixels

    corrmat =corr(data,data);    
    %Get the average correlation for each lag
    for i = 1:length(unique(lagmat))
        metric(k,i)=nanmean(corrmat(lagmat==i));        
    end
    
    %calculate the half life of the decay just using the active timepoints
    %only calculate tau for motifs that have at least 12 timepoints,
    %otherwise it's a bad fit for the eponential funciton
    temp = metric(k,:);
    temp(isnan(temp))=[];
    temp = temp(2:end); %ignore autocorelation
    if length(temp)==12
        x = (1:75:(numel(temp))*75)';
        y = temp';
        f = fit(x,y,'exp1','StartPoint',[0,0]);
        tau_hl(k) = -1*(log(2)/f.b);
    end
end
%Remove autocorrelations
metric(:,1)=[];

%report the number
fprintf('\nYou removed %d tau for incompleteness',sum(isnan(tau_hl)));
fprintf('\nYou removed  %d metric for no active timepoints',sum(all(isnan(metric),2))); 

rmv_tau = sum(isnan(tau_hl));
rmv_metric = sum(all(isnan(metric),2));
%Remove rows that are all NaN
metric(all(isnan(metric),2),:)=[];
tau_hl(isnan(tau_hl))=[];

%remove tau that are all nan

end %function












