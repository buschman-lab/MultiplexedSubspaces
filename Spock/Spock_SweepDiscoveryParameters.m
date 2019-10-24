function Spock_SweepDiscoveryParameters(block)

%Set path for spock to all repository (one step above cur dir)
cd ..
setpath(genpath(pwd));

save_dir = '/TrainRepitoires/TrainingFit_WindowLengthSweep';
save_file_name = sprintf('block_%d',block);

%% Body
%Parameter Combinations
% L = [1:12:65];

%Perform CNMF sweep
paramsweep = [];
for cur_val = [1:12:65]
    [paramsweep(cur_k).ExpVar_frame,paramsweep(cur_k).ExpVar_all,~,~,~,~,~,~,paramsweep(cur_k).numFactors,~,~,opts] = ...
        Discover_Motifs(...
            'block',block,...
            'K',50,...
            'L',cur_val);
end

filename = [opts.bucket opts.base save_dir save_file_name '.mat'];
save(filename,'paramsweep','opts','-v7.3');

end
