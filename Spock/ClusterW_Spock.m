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
nanpxs_all = cellfun(@(x) x.nanpxs, temp,'UniformOutput',0);
gp = loadobj(feval(parameter_class));  

%If different animals then take shared pixels
if numel(unique(cellfun(@(x) numel(x),nanpxs_all,'UniformOutput',1)))>1
   [W, nanpxs] = CrossAnimalMask(W_orig,nanpxs_all);
else
    W = W_orig;
    nanpxs = nanpxs_all{1}; 
end
fprintf('\n\t Generating Basis Motifs');
[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs,shift] = ClusterW(W,gp,nanpxs);

save([save_dir, filesep, file_name, '_basis_motifs.mat'],'W_basis', 'kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat','lag_mat','lags','nanpxs','gp','-v7.3')
saveCurFigs(handles,'-dpng',[file_name '_ClusteringMotifs'],save_dir,0); close all;

%create basis motifs per mouse
mouseid = MouseNumFromPath(file_list,'Mouse_'); unique_mice = unique(mouseid);
gp.clust_community_fraction = 0.75; gp.clust_removepad=1;
[W_mouse,core_comm_size,mouseid_all,cluster_idx_all] = WithinMouseBasisMotifs(W_orig,nanpxs_all,gp,mouseid,cluster_idx,idx_knn,lag_mat,lags,shift);
%Make sure that the timecourse of the basis motifs are alligned across mice
%(i.e. peak at same time)
if sum(diff(cellfun(@(x) size(x,2),W_mouse,'UniformOutput',1)))==0 %no mouse-specific noise motifs 
   W_mouse = AllignW_AcrossMice(W_mouse);
end
for i = 1:numel(W_mouse)
    W_basis = W_mouse{i};    
    save([save_dir, filesep, file_name, 'permouse_basis_motifs', num2str(unique_mice(i)), '.mat'],...
        'W_basis', 'mouseid','core_comm_size','kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat',...
        'lag_mat','lags','nanpxs','mouseid_all','cluster_idx_all','nanpxs_all','-v7.3')
end

fprintf('\n\t Done Saving Basis Motifs - Saved');
end