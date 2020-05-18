function FitMotifs_Spock(fn,save_fn,chunk,varargin)
%Camden MacDowell - timeless
%Spock shell for calling the FitMotifs function and saving off the results.

if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
end

%get figure names 
[fig_path, fig_name] = fileparts(save_fn);

%load the training and test data
temp = load(fn,'data_train','data_test','nanpxs');
nanpxs = temp.nanpxs; %store to save off for easier access later
data_train = squeeze(temp.data_train(:,chunk,:));
data_test = squeeze(temp.data_test(:,chunk,:));

[stats_test, stats_train, w, h, gp, fh] = FitMotifs(data_test,data_train);

saveCurFigs(fh,'-dpng',[fig_name '_lambdafit'],fig_path,0); close all;

%analyze the residuals of the training
handles = AnalyzeResiduals(data_train,helper.reconstruct(w,h),nanpxs,[]);
saveCurFigs(handles,'-dpng',[fig_name '_Res_train'],fig_path,0); close all;

%save off the data in the scratch directory and the nanpxs
fprintf('\n\tSaving data')

%save off the information in the scratch directory
save(save_fn,'stats_test','stats_train','w','h','nanpxs','gp','-v7.3')

fprintf('\n\tDONE')

end
















