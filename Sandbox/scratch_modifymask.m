

data = conditionDffMat(data_norm(:,1)',nanpxs);

[coef,score] = pca(data_norm(:,1:10000));
score = conditionDffMat(score',nanpxs);

imagesc(score(:,:,1)); axis square
manual_mask = impoly;
wait(manual_mask);               
manual_mask = manual_mask.createMask();
mask = manual_mask; 

data = reshape(data,[68*68,size(data,3)]);
data(mask(:),:) = NaN;
data = reshape(data,68,68,size(data,2));
[temp_data, temp_nanpxs] = conditionDffMat(new_mask);

[coef,score] = pca(temp_data(1:10000,:)');
score = conditionDffMat(score',temp_nanpxs);
close all; figure; hold on; for i = 1:20; subplot(4,5,i); imagesc(score(:,:,i)); axis square; end

%build a mask from the difference
orig_mask = conditionDffMat(nanmean(data_norm,2)',nanpxs);
new_mask = conditionDffMat(nanmean(temp_data,1),temp_nanpxs);
new_mask(49,9) = NaN;
new_mask(46:47,8) = nanmean(new_mask(:));
new_mask(47:49,9) = nanmean(new_mask(:));
new_mask(49,10) = nanmean(new_mask(:));
[temp_mask, temp_nanpxs] = conditionDffMat(repmat(new_mask,1,1,3));
new_mask = conditionDffMat(temp_mask(1,:),temp_nanpxs);

new_mask(~isnan(new_mask))=1;
new_mask(isnan(new_mask))=0;
orig_mask(~isnan(orig_mask))=1;
orig_mask(isnan(orig_mask))=0;
mask = orig_mask-new_mask;
close all; imagesc(mask)
mask(49,38)=1;
mask([58,58,56,56,56,57,57],[19,21,44,45,48,25,26])=1;
mask(59,17:28)=1;
mask(58,28)=1;


mask(58,29)=1;
mask(50,62)=1;
mask(9,28)=1;

load('Mouse331_RestingState_NP_06_11_2021_1dff_combined_processed','data_test','nanpxs');

figure; hold on; 
for i = 1:numel(file_list_preprocessed)
    load(file_list_preprocessed{i},'data_test','nanpxs');
    data = conditionDffMat(squeeze(data_test(:,:,1))',nanpxs);
    data = reshape(data,[68*68 size(data,3)]);
    data(mask(:)==1,:) = NaN;
    data = reshape(data,[68 68 size(data,2)]);
    subplot(3,2,i); imagesc(nanmean(data,3)); axis square
%     data = conditionDffMat(squeeze(data_test(:,:,1))',nanpxs);
    
end
    

