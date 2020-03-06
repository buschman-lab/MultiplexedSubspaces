temp = load('ClusteredDPs_rawdata.mat','W_all');
w = AlignMotifs(temp.W_all);

%get weight of each motif. 
amp_weight_norm = NaN(size(w,2),size(w,3));
for i = 1:size(w,2)
    temp = squeeze(w(:,i,:));
    temp(temp<=eps)=NaN;
    %intensity
    amp = nanmean(temp);
    num_active = sum(~isnan(temp),1);
    amp_weight_norm(i,:) = (amp.*num_active);%/max(amp.*num_active);
end

%plot the average intesn
figure; hold on;
x = (1:1:size(w,2))*75;
sem(amp_weight_norm');
shadedErrorBar(x,nanmean(amp_weight_norm)

%get the time to peak and the peak width (at halfway for each motif)





% 
% amp_weight_norm = NaN(size(w,2),size(w,3));
% for i = 1:size(w,2)
%     temp = squeeze(w(:,i,:));
%     temp(temp<=eps)=NaN;
%     %intensity
%     amp = nanmean(temp);
%     num_active = sum(~isnan(temp),1);
%     amp_weight_norm(i,:) = (amp.*num_active);%/max(amp.*num_active);
% end

% weight = NaN(size(w,2),size(w,3));
% for i = 1:size(w,2)
%     temp = squeeze(w(:,i,:));
%     %intensity
%     weight(i,:) = mean(temp);
% end