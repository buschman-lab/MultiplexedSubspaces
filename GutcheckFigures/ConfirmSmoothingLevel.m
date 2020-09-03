function ConfirmSmoothingLevel(data_train,nanpxs,W_basis)
%Confirms that the smoothing of data_train approximately matches the
%smoothing inherent in the basis motifs (because they are averages of
%multiple motifs/animals). 

%data_train is px x time, non nanpxs
W = W_basis; 

data_train = reshape(conditionDffMat(data_train',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));    

%Only work on shared pixels
[~, bad_pxl] = SharedPixels(W,data_train);
data_train(bad_pxl,:) = [];
W(bad_pxl,:,:) = [];

temp = {};
for i = 1:size(W,2)
    temp{i} = nanvar(squeeze(W(:,i,:)),[],1);
end
temp = [temp{:}];


%Loop through different smoothing levels
sm_sem = [];
sm_avg = [];
kernel = [0.1 0.25 0.5 0.75 1];
for i = 1:numel(kernel)
    smt_kernel = [kernel(i) kernel(i)];
    data_train_smooth = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_train_smooth(bad_pxl,:) = [];
    sm_sem(i) = sem(nanvar(data_train_smooth,[],1),2);
    sm_avg(i) = nanmean(nanvar(data_train_smooth,[],1),2);
end
    
y_err = [sem(temp,2),sem(nanvar(data_train,[],1),2),sm_sem];
y_mean = [nanmean(temp,2),nanmean(nanvar(data_train,[],1),2),sm_avg];

figure; hold on; 
barwitherr(y_err,y_mean)
set(gca,'xticklabel',[0,0,0,kernel])
ylabel('Spatial Variance');
xlabel('smooth kernel (starts with basis)')
title('determining the correct smoothing level'); 