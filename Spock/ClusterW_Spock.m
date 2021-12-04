function ClusterW_Spock(file_path,file_name, save_dir,parameter_class)
%Camden MacDowell - timeless
%spock shell to run motif clustering. 
% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

file_list = GrabFiles([file_name, '\w*chunk\w*.mat'],0,{file_path});
%load all the basis motifs
temp = cellfun(@(x) load(x,'w','nanpxs'),file_list,'UniformOutput',0);
W_orig = cellfun(@(x) x.w, temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.nanpxs, temp,'UniformOutput',0);
gp = loadobj(feval(parameter_class));  

%If different animals then take shared pixels
if numel(unique(cellfun(@(x) numel(x),nanpxs,'UniformOutput',1)))>1
   [W, nanpxs_all] = CrossAnimalMask(W_orig,nanpxs);
else
    W = W_orig;
    nanpxs_all = nanpxs{1}; 
end
fprintf('\n\t Generating Basis Motifs');
[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterW(W,gp,nanpxs_all);

save([save_dir, filesep, file_name, '_basis_motifs.mat'],'W_basis', 'kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat','lag_mat','lags','nanpxs','gp','-v7.3')
saveCurFigs(handles,'-dpng',[file_name '_ClusteringMotifs'],save_dir,0); close all;

%create basis motifs per mouse
mouseid = MouseNumFromPath(file_list,'Mouse_');
[W_mouse,core_comm_size] = WithinMouseBasisMotifs(W_orig,nanpxs,gp,mouseid,cluster_idx,idx_knn,lag_mat,lags);
save([save_dir, filesep, file_name, 'permouse_basis_motifs.mat'],'W_mouse', 'mouseid','core_comm_size','nanpxs','gp','-v7.3')

fprintf('\n\t Done Saving Basis Motifs - Saved');
end