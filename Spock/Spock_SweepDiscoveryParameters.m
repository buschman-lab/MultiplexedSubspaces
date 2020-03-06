function Spock_SweepDiscoveryParameters(block)

%Set path for spock to all repository
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

save_dir = '/TrainRepitoires/TrainingFit_WindowLengthSweep_FullData/';
save_file_name = sprintf('block_%d',block);


%% Body
%Parameter Combinations
% L = [1:12:65];

%Perform CNMF sweep
paramsweep = [];
COUNT = 0 ; 
for cur_val = [1:3:13,26:13:65]
    fprintf('\n\t Working on L = %d\n',cur_val);
    COUNT = COUNT+1;
    [paramsweep(COUNT).ExpVar_frame,paramsweep(COUNT).ExpVar_all,...
        ~,~,paramsweep(COUNT).H,~,~,~,paramsweep(COUNT).numFactors,...
        paramsweep(COUNT).W,paramsweep(COUNT).w,opts] = ...
        Discover_Motifs(...
            'block',block,...
            'K',50,...
            'L',cur_val);
     paramsweep(COUNT).L = cur_val; 
     paramsweep(COUNT).K = 50; 
end

save_dir = [opts.bucket opts.base save_dir];
if ~exist(save_dir)
    mkdir(save_dir);
end


filename = [save_dir save_file_name '.mat'];
save(filename,'paramsweep','opts','-v7.3');

end
