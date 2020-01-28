function Spock_ClusterBehavior(block)

addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
load('data_pheno');
k = [50,30,60,90,120,500];
[ovr_comm, ovr_Q] = PhenographSimple(features_downsampled, 'k', k(block));
save(sprintf('block%d',block'),'ovr_comm','ovr_Q','k','-v7.3');

end