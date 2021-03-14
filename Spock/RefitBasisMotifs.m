function RefitBasisMotifs(fn, basis_dir,chunk,parameter_class,save_dir)

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
end

gp = loadobj(feval(parameter_class)); 

%load the basis motifs
temp = load([basis_dir filesep 'basismotifs.mat'],'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W = temp.W_basis;
W(:,temp.noise_clusters,:)=[];

%load the data_train and the data_test data
fprintf('\n\tLoadings Data\n'); 
temp = load(fn,'data_test','data_train','nanpxs');
data_train = temp.data_train(:,:,chunk);
data_test = temp.data_test(:,:,chunk);
nanpxs = temp.nanpxs; 
[~,name] = fileparts(fn);
name = [name, sprintf('refitchunk_%d',chunk)];

%recondition, smooth, and flatten
if numel(gp.smt_kernel)>2 %autodetermine appropriate smoothing value
    fprintf('\n\t Autofitting Smoothing Value');
    gp.smt_kernel = AutoFitSmoothingLevel(cat(2,data_train,data_test),nanpxs,W,gp);
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
elseif numel(gp.smt_kernel)==2
    fprintf('\n\t Set Smoothing Value');
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
else
    fprintf('\n\t No Smoothing');
    data_train = reshape(conditionDffMat(data_train',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(conditionDffMat(data_test',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
end

%Only work on shared pixels
[~, bad_pxl] = SharedPixels(W,cat(2,data_test,data_train));
data_train(bad_pxl,:) = [];
data_test(bad_pxl,:) = [];
W(bad_pxl,:,:) = []; 
    
%fit motif weightings to data
[w,H] = fpCNMF(data_train,'non_penalized_iter',...
    gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
    'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
    'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
stats_refit = CNMF_Stats(w,H,data_train,0);
stats_refit.smoothingkernel = gp.smt_kernel;

%get residuals
residuals = data_train-tensor_convolve(w,H);

%save off
save([save_dir filesep name 'train.mat'],'w','H','stats_refit','bad_pxl','residuals');

[w,H] = fpCNMF(data_test,'non_penalized_iter',...
    gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
    'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
    'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
stats_refit = CNMF_Stats(w,H,data_test,0);
stats_refit.smoothingkernel = gp.smt_kernel;

%get residuals
residuals = data_test-tensor_convolve(w,H);

save([save_dir filesep name 'test.mat'],'w','H','stats_refit','bad_pxl','residuals');


%save off
