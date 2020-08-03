%VPA_Pipeline
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%Add paths
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis'));

%Load data location
temp = load('AllOriginalDataFileList.mat');
file_list = temp.file_list;

%save location 
save_dir ='Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\DeconvolutionParameterSweep';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end


%% Post processing and fit motifs

job_id = cell(1,numel(file_list));
rng('default');
idx = randperm(numel(file_list),5); %grab 5 random files
for cur_file = 1:numel(idx)     
    %deconvolve and split the data
    script_name = WriteBashScript(sprintf('%d',idx(cur_file)),'SweepDeconvolutionParameters',{idx(cur_file),...
        ConvertToBucketPath(save_dir)},...
        {'%d',"'%s'"},...
        'sbatch_time',1200,'sbatch_memory',12,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/GutcheckFigures/");
    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory of the bash script
        sprintf('sbatch %s',script_name)]);         
    %get job id
    job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
end

ssh2_close(s_conn);
clear username password sconn


