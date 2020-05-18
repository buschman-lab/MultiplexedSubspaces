function FitBasisMotifs_Spock(file_path,file_name, save_dir)
%Camden MacDowell - timeless
%spock shell to run motif clustering. 

% Add repository paths (you stark in dynamic scripts)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

file_list = GrabFiles([file_name, '\w*chunk\w*.mat'],0,{file_path});

gp=general_params; 
fprintf('\n\t Generating Basis Motifs');
[basis_motifs, kval, ovr_q, idx_louvain, idx_knn, tcorr_mat, handles] = ClusterMotifs(file_list,gp);
save([save_dir, filesep, file_name, '_basis_motifs.mat'],'basis_motifs', 'kval', 'ovr_q', 'idx_louvain', 'idx_knn', 'tcorr_mat','gp','-v7.3')
saveCurFigs(handles,'-dpng',[file_name '_ClusteringMotifs'],save_dir,0); close all;

fprintf('\n\t Done Saving Basis Motifs - Saved');
%% Now add the refitting code




end
