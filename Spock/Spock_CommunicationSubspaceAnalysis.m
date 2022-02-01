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

cur_rec=1;
%loop through each motif
for cur_m = 1:15
    script_name = WriteBashScript(parameter_class,sprintf('%d',1),'CommunicationSubspace_MotifTriggered',{cur_m,cur_rec},...
        {'%d','%d'},...
        'sbatch_time',480,'sbatch_memory',24,...
        'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/CommunicationSubspace/");

    ssh2_command(s_conn,...
    ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
    sprintf('sbatch %s',script_name)]);  
end


%close out connection
ssh2_close(s_conn);
clear username password s_conn

end %function end
