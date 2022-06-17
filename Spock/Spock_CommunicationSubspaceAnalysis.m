function Spock_CommunicationSubspaceAnalysis()
%Camden MacDowell
parameter_class = 'general_params_corticaldynamics';
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
end

% cur_rec=1;
muaflag = 1; 
% ccaflag=3; %1 cca local 2 cca full, 0 is rrr local and FA, 3 is rrr full, 4 is rrr_full reversed
%loop through each motif
for cur_rec = 1:6
    for cur_m = 1:15
        for ccaflag = 5:6 %[0,2,3,4]
        if ccaflag==0
            t=1080;
        else
            t=480;
        end        
%         script_name = WriteBashScript(parameter_class,sprintf('%d',1),'CommunicationSubspace_MotifTriggered_meansubtract',{cur_m,cur_rec,muaflag,ccaflag},...
        script_name = WriteBashScript(parameter_class,sprintf('%d',1),'CommunicationSubspace_MotifTriggered',{cur_m,cur_rec,muaflag,ccaflag},...
            {'%d','%d','%d','%d'},...
            'sbatch_time',t,'sbatch_memory',24,...
            'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/CommunicationSubspace/");

        ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);  
        end
    end
end


%% to run permutations on the ridge
%need to manually run each motif, otherwise too many jobs and it hangs
for cur_rec = 1:6
    for cur_m = 14 
        script_name = WriteBashScript(parameter_class,sprintf('%d',1),'RidgeRegressionPermutations',{'$SLURM_ARRAY_TASK_ID',cur_rec,cur_m},...
            {'%s','%d','%d'},...
            'sbatch_time',5,'sbatch_memory',12,'sbatch_array','1-1000',...
            'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/CommunicationSubspace/");

        ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);            
    end
end

%% to run LOO analysis
%load all the files
[fn,~] = GrabFiles('\w*regRRR_muaflag1_GROUPEDmotif\d*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace'});
for i = 1:numel(fn)
    script_name = WriteBashScript(parameter_class,sprintf('%d',1),'LeaveOneOutSubspaceFit',{ConvertToBucketPath(fn{i})},...
        {"'%s'"},...
        'sbatch_time',60,'sbatch_memory',12,'sbatch_array','1',...
        'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/CommunicationSubspace/");

    ssh2_command(s_conn,...
    ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
    sprintf('sbatch %s',script_name)]);      
end

%close out connection
ssh2_close(s_conn);
clear username password s_conn


end %function end
