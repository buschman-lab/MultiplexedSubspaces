function Spock_ClusterBehavior(block,k)

load('data_pheno');
[ovr_comm, ovr_Q] = PhenographSimple(features_downsampled, 'knn', k);
save(sprintf('block%d',block'),'ovr_comm','ovr_Q');

end