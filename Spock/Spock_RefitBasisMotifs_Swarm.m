function Spock_RefitBasisMotifs_Swarm(file_list,basis_dir,job_id,s_conn,parameter_class,save_dir)
%Camden MacDowell

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));    
end

%generate swarm
for cur_file = 1:numel(file_list)     
    temp = load(file_list{cur_file},'data_test');    
    for cur_chunk = 1:size(temp.data_test,3)
        script_name = WriteBashScript(sprintf('%d',1),'RefitBasisMotifs',{ConvertToBucketPath(file_list{cur_file}),...
            ConvertToBucketPath(basis_dir),cur_chunk,parameter_class,ConvertToBucketPath(save_dir)},...
        {"'%s'","'%s'",'%d',"'%s'"},...
        'sbatch_time',59,'sbatch_memory',10,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

        if ~isempty(job_id)
            ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch --dependency=afterok:%s %s',job_id,script_name)]); 
        else
            ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]);  
        end
    end
end

