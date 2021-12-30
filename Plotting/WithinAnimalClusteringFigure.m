function WithinAnimalClusteringFigure()
%Basis motif figure/within animal clustering
%Camden MacDowell - timeless

cd('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\MotifDiscovery');
basis = 'Mou_basis_motifs.mat';
mouse_fn = {'Moupermouse_basis_motifs331.mat','Moupermouse_basis_motifs332.mat','Moupermouse_basis_motifs334.mat'};

load(basis); 

%plot the clustering matrix
