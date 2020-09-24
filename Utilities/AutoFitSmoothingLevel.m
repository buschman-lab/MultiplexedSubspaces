function kernel = AutoFitSmoothingLevel(data,nanpxs,W,gp)
%Camden MacDowell - timeless
%Matches the smoothing level of data to the spatial variance of W_basis.
%This is necessary because the basis motifs are averages of multiple motifs
%across animals/sessions so spatial variance is smoothed

%data is px x time, non nanpxs
data = reshape(conditionDffMat(data',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data,2));    

%Only work on shared pixels
[~, bad_pxl] = SharedPixels(W,data);
data(bad_pxl,:) = [];
W(bad_pxl,:,:) = [];

%Get the spatial variance of W. Use a loop so no weird reshape mistakes
w_var = cell(1,size(W,2));
for i = 1:size(W,2)
    w_var{i} = nanvar(squeeze(W(:,i,:)),[],1);
end
w_var = [w_var{:}];

%Loop through different smoothing levels
sm_avg = zeros(1,numel(gp.smt_kernel));
for i = 1:numel(gp.smt_kernel)
    smt_kernel = [gp.smt_kernel(i) gp.smt_kernel(i)];
    data_train_smooth = reshape(SpatialGaussian(conditionDffMat(data',find(bad_pxl==1)),smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data,2));
    data_train_smooth(bad_pxl,:) = [];
    sm_avg(i) = nanmean(nanvar(data_train_smooth,[],1),2);  
end

%get the closest value
[~,idx] = min(abs(sm_avg-nanmean(w_var,2)));

kernel = [gp.smt_kernel(idx),gp.smt_kernel(idx)];

end %function end

% % %% figure; 
% y_mean = [nanmean(w_var,2),nanmean(nanvar(data,[],1),2),sm_avg];
% 
% figure; hold on; 
% plot(y_mean)
% ylabel('Spatial Variance');
% xlabel('smooth kernel (starts with basis)')
% title('determining the correct smoothing level')




