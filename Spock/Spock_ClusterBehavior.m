function Spock_ClusterBehavior(block,k)

addpath(getpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
load('data_pheno');
[ovr_comm, ovr_Q] = PhenographSimple(features_downsampled, 'knn', k);
save(sprintf('block%d',block'),'ovr_comm','ovr_Q');

end