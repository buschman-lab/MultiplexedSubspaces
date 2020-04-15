%PREPROCESSING PIPELINE
% Camden MacDowell 2020
% User selects folders of recordings for one animal and does manual steps
% Then automatically generates spock bash scripts to run and
% send a dependency to spock to combined all dffs in each folder
% in order. 

%% Manual steps
save_dir = 'Z:\Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Widefield_Preprocessed\'; %Final target save directory for the preprocesssed and combined files
if ~exist(save_dir)
    mkdir(save_dir);
end

%select folders to process and grab the first file from each rec.
[file_list,folder_list] = GrabFiles('Pos0.ome.tif');

%% Manual Configure Preprocessing
%configure preprocessing options
opts = ConfigurePreProcessing('crop_w',300,'vasc_std',2,'save_uncorrected',0);
 
%Grab reference images for each. Preload so no delay between loop.
ref_imgs = cellfun(@(x) GetReferenceImage(x,opts.fixed_image),...
    file_list, 'UniformOutput',0);

%loop through reference images, register, and apply manual changes
for cur_fold = 1:numel(folder_list)
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
        [~,fn] = fileparts(folder_list{cur_fold});
        saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),[save_dir 'PreprocessingFigures'],0); %close all;       
        close all;
    end
    %save off the options to each folder
    save([folder_list{cur_fold} filesep 'prepro_log'],'prepro_log')
end


%% Spock Preprocessing 
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%gather recordings 
for cur_fold = 1:numel(folder_list)
    [file_list,~] = GrabFiles('.tif',0,folder_list(cur_fold));
    [opts_list,~] = GrabFiles('prepro_log.mat',0,folder_list(cur_fold)); 
    
    %Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list));
    input_type = {"'%s'","'%s'"};
    for cur_file = 1:numel(file_list)
        input_val = {ConvertToBucketPath(file_list{cur_file}), ConvertToBucketPath(opts_list{1})};
        script_name = WriteBashScript(sprintf('%d_%d',cur_fold,cur_file),'Spock_Preprocessing_Pipeline',input_val,input_type);

        %Run job
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]); 
        
        %get job id
        job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
        if cur_file ~=numel(file_list)
            job_id{cur_file} = [job_id{cur_file} ','];
        end
    end    
       
    %Once each folder is done, combine all the dffs in that folder
    input_type = {"'%s'","'%s'"};
    input_val = {ConvertToBucketPath(folder_list{cur_fold}), ConvertToBucketPath(save_dir)};    
    script_name = WriteBashScript(sprintf('%d_%d_combine',cur_fold,cur_file),'Spock_CombineStacks',input_val,input_type);    
    
    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]);
    
end



