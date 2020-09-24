function ClusterW_Spock(file_path,file_name, save_dir)
%Camden MacDowell - timeless
%spock shell to run motif clustering. 
% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

file_list = GrabFiles([file_name, '\w*chunk\w*.mat'],0,{file_path});
%load all the basis motifs
temp = cellfun(@(x) load(x,'w','nanpxs'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.w, temp,'UniformOutput',0);
nanpxs = cellfun(@(x) x.nanpxs, temp,'UniformOutput',0);
gp=general_params; 
fprintf('\n\t Generating Basis Motifs');
[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterW(W,gp,nanpxs{1});

save([save_dir, filesep, file_name, '_basis_motifs.mat'],'W_basis', 'kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat','lag_mat','lags','nanpxs','gp','-v7.3')
saveCurFigs(handles,'-dpng',[file_name '_ClusteringMotifs'],save_dir,0); close all;

fprintf('\n\t Done Saving Basis Motifs - Saved');
end