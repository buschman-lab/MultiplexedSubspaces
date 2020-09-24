function FitMotifs_Spock(fn,save_fn,chunk,parameter_class)
%Camden MacDowell - timeless
%Spock shell for calling the FitMotifs function and saving off the results.

if nargin <3
    chunk = [];
end

if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

%load the training and test data
gp = loadobj(feval(parameter_class)); 
temp = load(fn,'data_train','data_test','nanpxs');
nanpxs = temp.nanpxs; %store to save off for easier access later

if ~isempty(chunk)
    data_train = squeeze(temp.data_train(:,:,chunk));
    data_test = squeeze(temp.data_test(:,:,chunk));
else
    data_train = temp.data_train;
    data_test = temp.data_test;
end

if chunk == 1 %make some example figures
    [w, h, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,gp,1);
    %get figure names 
    [fig_path, fig_name] = fileparts(save_fn);
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-dpng',[fig_name '_examplefit'],fig_path,0); close all;
    %analyze the residuals of the training
    handles = AnalyzeResiduals(data_train,tensor_convolve(w,h),nanpxs,[]);
    saveCurFigs(handles,'-dpng',[fig_name '_Res_train'],fig_path,0); close all;

else %no figures
    [w, h, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,gp,0);
end  

%save off the data in the scratch directory and the nanpxs
fprintf('\n\tSaving data')

%save off the information in the scratch directory
save(save_fn,'stats_test','stats_train','w','h','nanpxs','-v7.3')

fprintf('\n\tDONE')

end
















