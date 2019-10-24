% take the dff and perform PCA on it. 

dff = load('C:\Users\macdo\Documents\ExampleDataFromLab\Making Raw Movies\Mouse4757_RestingState_1\Mouse4757_RestingState_Processing\DFFCombined_filt_lin_bin_smooth.mat');
nanpxs = dff.nanpxs_array{1};
dff = dff.dff_norm_lin_array{1}; 

[time_coef, time_score, time_latent] = pca(dff);  %PCA in timepoints


[coef, score, latent] = pca(dff');  %PCA in space
temp = conditionDffMat(coef',nanpxs);



