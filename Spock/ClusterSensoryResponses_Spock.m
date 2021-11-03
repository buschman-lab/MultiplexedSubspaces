function ClusterSensoryResponses_Spock(file_path,file_name, save_dir,type)
%Camden MacDowell - timeless
%spock shell to run motif clustering. 
% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis/'));
    cd('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock');
end

file_list = GrabFiles(['Mouse9036','.*',num2str(type),'.mat'],0,{file_path}); %could add file_name here with 'vis1' or something.
%load all the basis motifs
temp = cellfun(@(x) load(x,'data','nanpxs'),file_list,'UniformOutput',0);
W = cellfun(@(x) x.data(12:12+12,:)', temp,'UniformOutput',0);

%remove pixels with zero variance across all
nanpxs = find(nanvar(cat(2,W{:}),[],2)<=eps);

%concatenate along the second dimension
W = cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),W,'UniformOutput',0);
W = cat(2,W{:});
W(nanpxs,:,:)=[];

gp=general_params_sensoryResponses; 
fprintf('\n\t Generating Sensory Response Clusters');
[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterSensoryResponses(W,gp,nanpxs);

save([save_dir, filesep, file_name, '_',num2str(gp.clust_knn),'_responseClusters.mat'],'W_basis', 'kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat','lag_mat','lags','nanpxs','gp','-v7.3')
saveCurFigs(handles,'-dpng',[file_name '_responseClusters',num2str(gp.clust_knn)],save_dir,0); close all;

fprintf('\n\t Done Saving Basis Motifs - Saved');
end