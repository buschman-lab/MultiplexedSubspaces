function Spock_STDecomposition(block,save_dir)

%Set path for spock to all repository (one step above cur dir)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

if nargin<2
    save_dir = '/TrainRepitoires/SpaceTimeNMF/';
end
save_file_name = sprintf('block_%d',block);

%% Body
[X_recon,err,Wt,A,Ws, fit_to_data, opts] = DecomposeMotifs('block',block);

%Save off
save_dir = [opts.bucket opts.base save_dir];
if ~exist(save_dir)
    mkdir(save_dir);
end

filename = [save_dir save_file_name '.mat'];
save(filename,'X_recon','err','Wt','A','Ws', 'fit_to_data','-v7.3');

fprintf('Done Saving Data to %s',filename');

end
