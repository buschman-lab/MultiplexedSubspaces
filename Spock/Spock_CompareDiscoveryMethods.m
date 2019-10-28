function Spock_CompareDiscoveryMethods(block,save_dir)

%Set path for spock to all repository (one step above cur dir)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

if nargin<2
    save_dir = '/TrainRepitoires/TrainingFit_CompareDiscoveryMethods/';
end
save_file_name = sprintf('block_%d',block);

%% Body

%Perform seqNMF
[seq.ExpVar_frame,seq.ExpVar_all,~,~,~,~,~,~,seq.numFactors,~,~,seq.opts] = ...
    Discover_Motifs('block',block,'K',50,'bucket','Z:/');

%Perform CNMF sweep
nmf = struct();
for cur_k = [1:1:25,50:25:300]
    [nmf(cur_k).ExpVar_frame,nmf(cur_k).ExpVar_all,~,~,~,~,~,~,nmf(cur_k).numFactors,~,~,nmf(cur_k).opts] = ...
        Discover_Motifs(...
            'block',block,...
            'K',cur_k,...
            'L',1,...
            'lambda',0,...
            'lambdaL1H',0,...
            'lambdaOrthoH',0,...
            'RemoveEmptyWs',0);
end

%Perform PCA 
[spca.ExpVar_all,spca.coeff,spca.score,spca.mu,spca.opts] = Discover_Networks_PCA('block',block);

%Save off
save_dir = [seq.opts.bucket seq.opts.base save_dir];
if ~exist(save_dir)
    mkdir(save_dir);
end

filename = [save_dir save_file_name '.mat'];
save(filename,'nmf','spca','seq','-v7.3');

end
