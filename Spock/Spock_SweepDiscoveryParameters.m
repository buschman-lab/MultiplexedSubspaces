function Spock_SweepDiscoveryParameters(block)

%Set path for spock to all repository (one step above cur dir)
cd ..
addpath(genpath(pwd));

save_dir = '/TrainRepitoires/TrainingFit_WindowLengthSweep/';
save_file_name = sprintf('block_%d',block);


%% Body
%Parameter Combinations
% L = [1:12:65];

%Perform CNMF sweep
paramsweep = [];
COUNT = 0 ; 
for cur_val = [1:12:65]
    COUNT = COUNT+1;
    [paramsweep(COUNT).ExpVar_frame,paramsweep(COUNT).ExpVar_all,~,~,~,~,~,~,paramsweep(COUNT).numFactors,~,~,opts] = ...
        Discover_Motifs(...
            'block',block,...
            'K',50,...
            'L',cur_val);
end

save_dir = [opts.bucket opts.base save_dir];
if ~exist(save_dir)
    mkdir(save_dir);
end


filename = [save_dir save_file_name '.mat'];
save(filename,'paramsweep','opts','-v7.3');

end
