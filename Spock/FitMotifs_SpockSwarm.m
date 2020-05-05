function [swarm_id,save_fn] = FitMotifs_SpockSwarm(fn,dependency_id,s_conn)
%Camden MacDowell
%once file fn is done (dependency), parallellizes the fitting of each data
%chunk in fn

gp = general_params;
    
%load the training data
temp = load(fn,'num_chunks');
num_chunks = temp.num_chunks/2; %since chunks is both training and testing

%generate swarm
swarm_id = cell(1,num_chunks);
save_fn = cell(1,num_chunks);
for i = 1:num_chunks
    [~, fn_temp] = fileparts(fn);  
    save_fn{i} = [gp.local_bucket gp.processing_intermediates fn_temp sprintf('_fit_chunk%d.mat',i)];
        
    script_name = WriteBashScript(sprintf('motifchunk%d',i),'FitMotifs',{ConvertToBucketPath(fn),ConvertToBucketPath(save_fn{i}),i},{"'%s'","'%s'","%d"},...
        'sbatch_time',500,'sbatch_memory',12,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Analysis_functions/");
    
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',dependency_id,script_name)]);
    
    swarm_id{i} = erase(response.command_result{1},'Submitted batch job ');
    if i ~=num_chunks
        swarm_id{i} = [swarm_id{i} ','];
    end
end