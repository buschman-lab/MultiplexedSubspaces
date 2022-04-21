function Spock_CorticalNetworks()
%Also works for 'behavioral networks';
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
% ccaflag=3; %1 cca local 2 cca full, 0 is rrr local and FA, 3 is rrr full, 4 is rrr_full reversed
%loop through each motif
for cur_rec = 5:6 %1:6
    for cur_m = 11 %[1,3:15]
        for cur_a = 1:8 %[1:8]
            for reverseFlag = 0:1
                script_name = WriteBashScript(parameter_class,sprintf('%d',1),'CorticalNetworks',{cur_rec,cur_m,cur_a,reverseFlag},...
                    {'%d','%d','%d','%d'},...
                    'sbatch_time',239,'sbatch_memory',24,...
                    'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/CommunicationSubspace/");

%                 script_name = WriteBashScript(parameter_class,sprintf('%d',1),'BehavioralNetworks',{cur_rec,cur_m,cur_a,reverseFlag},...
%                     {'%d','%d','%d','%d'},...
%                     'sbatch_time',239,'sbatch_memory',14,...
%                     'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/BehavioralAnalysisFunctions/");

                ssh2_command(s_conn,...
                ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
                sprintf('sbatch %s',script_name)]);  
            end
        end
    end
end

%close out connection
ssh2_close(s_conn);
clear username password s_conn


end %function end
