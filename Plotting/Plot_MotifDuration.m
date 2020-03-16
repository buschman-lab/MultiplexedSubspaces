function Plot_MotifDuration(fn)
% load data 
for i = 1:numel(fn)
    fprintf('\n\tloading_data %d',i);
   w = load(fn{i},'paramsweep');
end
%% plot with an error bar
figure; hold on;
for i = 1:size(fn)
    fprintf('\n\tloading_data %d',i);
    w = load(fn{i},'paramsweep');
%     w = data{i};
    amp_weight = NaN(size(w,2),size(w,3));
    for j = 1:size(w,2)
        temp = squeeze(w(:,j,:));
        temp(temp<=eps)=NaN;
        %intensity
        amp = nanmedian(temp);
        num_active = sum(~isnan(temp),1);
        amp_weight(j,:) = (amp.*num_active);%/max(amp.*num_active);
    end
    x = (1:1:size(w,3))*75;
    err = nanstd(amp_weight,[],1)./sqrt(sum(~isnan(amp_weight),1));
    shadedErrorBar(x,nanmedian(amp_weight),err)
end
title('Average Activity Across Motifs');
legend('4','7','10','13','26','39')

%plot without an errorbar
figure; hold on;
for i = 1:size(data)
    w = data{i};
    amp_weight = NaN(size(w,2),size(w,3));
    for j = 1:size(w,2)
        temp = squeeze(w(:,j,:));
        temp(temp<=eps)=NaN;
        %intensity
        amp = nanmedian(temp);
        num_active = sum(~isnan(temp),1);
        amp_weight(j,:) = (amp.*num_active);%/max(amp.*num_active);
    end
    x = (1:1:size(w,3))*75;
    plot(x,nanmedian(amp_weight))
end
title('Average Activity Across Motifs');
legend('4','7','10','13','26','39')

%plot without an errorbar
figure; hold on;

%% Plot a bar plot of the 




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