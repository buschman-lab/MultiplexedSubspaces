%PREPROCESSING PIPELINE
% Camden MacDowell 2020
% User selects folders of recordings for one animal and does manual steps
% Then automatically generates spock bash scripts to run and
% send a dependency to spock to combined all dffs in each folder
% in order. 
% Currently only supports running a single animal at once 

%% THINGS TO ADD: 
%ability to run from start to finish
%commented sections for rerun from each break point
%combine with analysis pipeline
%Future addition: Remove recordings flagged with artifacts

%% Manual steps
save_dir_preprocessed = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_Preprocessed\'; %Final target save directory for the preprocesssed and combined files
save_dir_preprocessed = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_Preprocessed\'; %Final target save directory for the preprocesssed and combined files
if ~exist(save_dir_preprocessed)
    mkdir(save_dir_preprocessed);
end

save_dir_motif_fits = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_MotifFits\'; %Target directory of basis motifs and temporal weightings
if ~exist(save_dir_motif_fits)
    mkdir(save_dir_motif_fits);
end

save_dir_figs = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\ProcessingFigures\'; %Target directory of basis motifs and temporal weightings
if ~exist(save_dir_figs)
mkdir(save_dir_figs);
end

% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%select folders to process and grab the first file from each rec.
[file_list_first_block,folder_list_raw] = GrabFiles('Pos0.ome.tif');

%configure preprocessing options
opts = ConfigurePreProcessing('crop_w',300,'vasc_std',2,'save_uncorrected',0);
%load general params (this is for anything after preprocessing)
gp = general_params;

%% Manual Processing

%Grab reference images for each. Preload so no delay between loop.
ref_imgs = cellfun(@(x) GetReferenceImage(x,opts.fixed_image),...
    file_list_first_block, 'UniformOutput',0);

%loop through reference images, register, and apply manual changes
for cur_fold = 1:numel(folder_list_raw)
    if cur_fold ==1
        %manual allignment 
        prepro_log = ManualAlignment(ref_imgs{cur_fold},opts);

        %mask vasculature and manual cleanup (optional)
        prepro_log = MaskVasculature(...
            prepro_log.cropped_alligned_img,prepro_log);
        
        close
        %no transformation
        prepro_log.tform = []; 
        prepro_log.output_size = [];
    else %new images
        %Register Reference Images to the first reference image          
        prepro_log = RegisterReferenceImages(ref_imgs{1},ref_imgs{cur_fold},prepro_log);               
        [~,fn] = fileparts(folder_list_raw{cur_fold});
        saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),save_dir_figs,0); %close all;       
        close all;
    end
    %save off the options to each folder
    save([folder_list_raw{cur_fold} filesep 'prepro_log'],'prepro_log')
end


%% Preprocessing. Results in a single hemo corrected, masked recording for each day
%gather recordings 
file_list_preprocessed = cell(1,numel(folder_list_raw));
for cur_fold = 1:numel(folder_list_raw)
    [file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold));
    [opts_list,~] = GrabFiles('prepro_log.mat',0,folder_list_raw(cur_fold)); 
    
    %Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    input_type = {"'%s'","'%s'"};
    for cur_file = 1:numel(file_list_raw)
        input_val = {ConvertToBucketPath(file_list_raw{cur_file}), ConvertToBucketPath(opts_list{1})};
        script_name = WriteBashScript(sprintf('%d_%d',cur_fold,cur_file),'Spock_Preprocessing_Pipeline',input_val,input_type,...
            'sbatch_time',15,'sbatch_memory',10);
        
        %Run job
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]); 
        
        %get job id
        job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
        if cur_file ~=numel(file_list_raw)
            job_id{cur_file} = [job_id{cur_file} ','];
        end
    end    
       
    %Once each folder is done, combine all the stacks and do hemocorrection
    [~,header] = fileparts(ConvertToBucketPath(folder_list_raw{cur_fold}));
    file_list_preprocessed{cur_fold} = [save_dir_preprocessed header 'dff_combined.mat']; 
    script_name = WriteBashScript(sprintf('%d_%d_combine',cur_fold,cur_file),'Spock_CombineStacks',{ConvertToBucketPath(folder_list_raw{cur_fold}),ConvertToBucketPath(file_list_preprocessed{cur_fold})},{"'%s'","'%s'"});    
    
    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]);    
end

%% Deconvolution/normalization and Motif Fitting. Results in cross validated motifs
% %optional restart: selected combined dff files
% file_list_preprocessed = GrabFiles('501\w*dff_combined.mat');
gp = general_params;

job_id = cell(1,numel(file_list_preprocessed));
file_list_motifs = cell(1,numel(file_list_preprocessed));
file_list_processed = cell(1,numel(file_list_preprocessed));
for cur_file = 1:numel(file_list_preprocessed)  
    [~, fn_temp] = fileparts(file_list_preprocessed{cur_file});      
    file_list_processed{cur_file} = [gp.local_bucket gp.processing_intermediates fn_temp '_processed.mat']; %path in temporary folder with the train/test split data
    
    %deconvolve and split the data    
    script_name = WriteBashScript(sprintf('%d',cur_file),'ProcessAndSplitData',{ConvertToBucketPath(file_list_preprocessed{cur_file}),...
        ConvertToBucketPath(file_list_processed{cur_file})},...
        {"'%s'","'%s'"},...
        'sbatch_time',5,'sbatch_memory',16,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);         
    %get job id
    temp_job_id = erase(response.command_result{1},'Submitted batch job ');
    
    %Fit motifs and cross-validate (parallellize by chunk)
    [swarm_id, swarm_motifs] = FitMotifs_SpockSwarm(file_list_processed{cur_file},temp_job_id,s_conn,gp);
    job_id{cur_file} = [swarm_id{:}]; 
    file_list_motifs{cur_file} = swarm_motifs;
end

%% Cluster motifs and refit to all data. Results in single set of motifs 
% file_list_motifs = GrabFiles('501\w*chunk\w*.mat');
% swarm_id = [];

[file_path_motifs, file_header_motifs] = fileparts(file_list_motifs{1});
header = file_header_motifs(1:regexp(file_header_motifs,'_','start','once'));

%Find Basis Motifs and Refit to the entire data set
script_name = WriteBashScript(sprintf('%d',1),'ClusterW_Spock',{ConvertToBucketPath(file_path_motifs),header,ConvertToBucketPath(save_dir_motif_fits)},...
    {"'%s'","'%s'","'%s'"},...
    'sbatch_time',600,'sbatch_memory',16,...
    'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

if ~isempty(swarm_id) % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[swarm_id{:}],script_name)]);    
else
   response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]); 
end


%% make figures, delete temporary data

%Run clean up script and close out connection
ssh2_close(s_conn);
clear username password sconn
%%



