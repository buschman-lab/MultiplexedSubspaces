%Analysis Pipeline
%Camden MacDowell 2020
% User selects preprocessed dffs. This loads, filters, and fits 

save_dir = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_MotifFits\'; %Final target save directory for the preprocesssed and combined files
if ~exist(save_dir)
    mkdir(save_dir);
end

%Todo
%autosplit by mice so that you just have to press go, and it runs

%select files to process
file_list = GrabFiles('dff_combined.mat');

%gut check: confirm is all recs are from the same animal
assert(unique(MouseNumFromFilename(file_list))~=1,'Error: user selected recordings from multiple mice'); 

%% Spock Motif Analysis
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%load general params (this is for anything after preprocessing)
gp = general_params;

%Filter, normalize, and break into testing and training chunks 
job_id = cell(1,numel(file_list));
for cur_file = 1:numel(file_list)  
    [~, fn_temp] = fileparts(file_list{cur_file});  
    save_split = [gp.local_bucket gp.processing_intermediates fn_temp '_processed.mat'];
    
    %process and split the data
    script_name = WriteBashScript(sprintf('%d',cur_file),'ProcessAndSplitData',{ConvertToBucketPath(file_list{cur_file}),ConvertToBucketPath(save_split)},{"'%s'","'%s'"},...
        'sbatch_time',5,'sbatch_memory',16,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);         
    %get job id
    temp_job_id = erase(response.command_result{1},'Submitted batch job ');
    
    %Fit motifs and cross-validate (parallellize by chunk)
    [swarm_id, save_motifs] = FitMotifs_SpockSwarm(save_split,temp_job_id,s_conn);
    job_id{cur_file} = [swarm_id{:}]; 
end


ssh2_close(s_conn);
clear username password sconn
%Cluster, and fit to entire recording (maybe add a lambda sweep?)




%












