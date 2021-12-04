function CompareMotifDiscoveryParamaters(data_dir)
%Camden MacDowell - timeless
%Performs motif discovery in all recordings in data dir using different
%parameters listed as temp in the ParameterClasses folder

if nargin<1; data_dir = {'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging'}; end 

[file_list_processed,~] = GrabFiles('\w*_processed.mat',0,data_dir); %select the preprocessed data (not the '_processed');


% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);


%load general params (this is for anything after preprocessing)
[parameter_class_all,~] = GrabFiles('temp\w*.m',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis\ParameterClasses'});
[~,parameter_class_all] = cellfun(@(x) fileparts(x),parameter_class_all,'UniformOutput',0);

for cur_p = 1:numel(parameter_class_all) 
    parameter_class = parameter_class_all{cur_p};
    gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string.
    file_list_motifs = cell(1,numel(file_list_processed));
    for cur_file = 1:numel(file_list_processed)      
        %Fit motifs and cross-validate (parallellize by chunk)
        [swarm_id, swarm_motifs] = FitMotifs_SpockSwarm(file_list_processed{cur_file},[],s_conn,parameter_class);
        job_id{cur_file} = [swarm_id{:}]; 
        file_list_motifs{cur_file} = swarm_motifs;
    end
    % Cluster Motifs
    save_dir_motif_fits = ['Z:\' gp.processing_intermediates];
    [file_path_motifs, ~] = fileparts(file_list_motifs{1}(1));
    header = ''; %leave blank to do by all animals %file_header_motifs(1:regexp(file_header_motifs,'_','start','once'));%to do by recording

    %Find Basis Motifs and Refit to the entire data set
    script_name = WriteBashScript(parameter_class,sprintf('%d',1),'ClusterW_Spock',{ConvertToBucketPath(file_path_motifs),header,ConvertToBucketPath(save_dir_motif_fits),parameter_class},...
        {"'%s'","'%s'","'%s'","'%s'"},...
        'sbatch_time',2879,'sbatch_memory',64,... %2879
        'sbatch_path',"/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/Spock/");

    if ~isempty(swarm_id) % Run job with dependency
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]);    
    else
       response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Projects/Cortical\ Dynamics/Cortical\ Neuropixel\ Widefield\ Dynamics/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]); 
    end
end

%Run clean up script and close out connection
ssh2_close(s_conn);
clear username password sconn

