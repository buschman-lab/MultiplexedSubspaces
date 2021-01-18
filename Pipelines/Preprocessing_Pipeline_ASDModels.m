% This pipeline walks the user through entire preprocessing, motif
% discovery, and motif fitting of mesoscale widefield calcium data. 
%It requires a computing cluster (e.g. spock at princeton). 
%Most of the complexity in this pipeline is to enable massive batch 
%processing of many animals and recordings at once. However the user can/should 
%Remove this complexity and get a better feel for the exact operations
%being performed by looking at the underlying code. 

%reach out to camdenm@princeton.edu with questions. 

%Core Code:
%%PreProcess.m: This move through each image in a stack, masks, spatially
%bins, and registers it to the reference image for that animal. 
%Also will register multiple sessions in the same animal 
%(if multiple sessions were selected). 
%%HemodynamicCorrection.m: separates the recording by wavelengths and
%corrects
%%makeDFF: Name says it all. 
%%ProcessAndSplitData: Splits into chunks and deconvolves data

%FINALLY: BE SURE TO CLOSE ALL SSH CONNECTIONS WHEN DONE!!!! OTHERWISE PNIHELP 
%GETS MAD. I STRONGLY RECCOMEND ADDING A WEEKLY REMINDER TO YOUR 
%CALENDAR TO KILL ALL OPEN CONNECTIONS WITH 'killall -u username'
%If you run this straight through it closes the connection at the very end 

% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%Add paths
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));

%configure preprocessing options
opts = ConfigurePreProcessing('crop_w',540,'vasc_std',2,'save_uncorrected',1,'fixed_image','first','spatial_bin_factor',4,'method_window',15);

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_asdmodels';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

%set up save directories
save_dir_processed = 'Z:\Projects\Cortical Dynamics\Mouse Models of Autism\SHANK Cohort 2020\Processed\'; %target savedirector
if ~exist(save_dir_processed,'dir')
    mkdir(save_dir_processed);
end


%% Manual Portion

%select folders to process and grab the first file from each rec.
%EXAMPLE DATA: Select 'folders' and then select 'Z:\Rodent Data\WideField Microscopy\ExampleData\Mouse431_10_17_2019\431-10-17-2019_1'
[file_list_first_stack,folder_list_raw] = GrabFiles('Pos0.ome.tif',1,{'Z:\Rodent Data\Wide Field Microscopy\ASD Models_Widefield'});

%Grab reference images for each. Preload so no delay between loop.
ref_imgs = cellfun(@(x) GetReferenceImage(x,opts.fixed_image),...
    file_list_first_stack, 'UniformOutput',0);

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
        %Register Reference Images to the first reference image try to do it automatedly, but backs up with manual allignment if necessary        
        prepro_log = RegisterReferenceImages(ref_imgs{1},ref_imgs{cur_fold},prepro_log);               
        [~,fn] = fileparts(folder_list_raw{cur_fold});
        saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),save_dir_processed,0); %close all;       
        close all;
    end
    %save off the options to each folder
    save([folder_list_raw{cur_fold} filesep 'prepro_log'],'prepro_log')
end


%% Preprocessing. Results in a single hemo corrected, masked recording for each day in the 'preprocessed' folder
%gather recordings
file_list_preprocessed = cell(1,numel(folder_list_raw));
for cur_fold = 1:numel(folder_list_raw)
    [file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold));
    [opts_list,~] = GrabFiles('prepro_log.m',0,folder_list_raw(cur_fold)); 
    
    %Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    for cur_file = 1:numel(file_list_raw)
        input_val = {ConvertToBucketPath(file_list_raw{cur_file}), ConvertToBucketPath(opts_list{1})};
        script_name = WriteBashScript(sprintf('%d_%d',cur_fold,cur_file),'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"},...
            'sbatch_time',60,'sbatch_memory',16);  %
        
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
    file_list_preprocessed{cur_fold} = [save_dir_processed header 'dff_combined.mat']; 
    script_name = WriteBashScript(sprintf('%d_combine',cur_fold),'Spock_CombineStacks',{ConvertToBucketPath(folder_list_raw{cur_fold}),ConvertToBucketPath(file_list_preprocessed{cur_fold}),parameter_class},{"'%s'","'%s'","'%s'"},...
        'sbatch_time',300,'sbatch_memory',128);    
    
    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]); 
    
%     % Run job with no dependency
%     response = ssh2_command(s_conn,...
%         ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
%         sprintf('sbatch %s',script_name)]); 
end

%if you want to see what the preprocessed data looks like then run
% InspectPreprocessedData(PreprocessedDataFilepath,'preprocessed')

%% Deconvolution, chunking, and Motif Fitting. Results in cross validated motifs in the MotifFits folder and Deconvolved and chunked data in the Preprocessed
%if you need to restart from this point: 
% [file_list_preprocessed,~] = GrabFiles('\w*combined.mat'); %select the preprocessed data (not the '_processed');

% NEED TO WAIT FOR ABOVE TO COMPLETE
job_id = cell(1,numel(file_list_preprocessed));
file_list_motifs = cell(1,numel(file_list_preprocessed));
file_list_processed = cell(1,numel(file_list_preprocessed));
for cur_file = 1:numel(file_list_preprocessed)  
    [~, fn_temp] = fileparts(file_list_preprocessed{cur_file});      
    file_list_processed{cur_file} = [save_dir_processed fn_temp '_processed.mat']; %path in temporary folder with the train/test split data
    
    %deconvolve and split the data    
    script_name = WriteBashScript(sprintf('%d',cur_file),'ProcessAndSplitData',{ConvertToBucketPath(file_list_preprocessed{cur_file}),...
        ConvertToBucketPath(file_list_processed{cur_file}),parameter_class},...
        {"'%s'","'%s'","'%s'"},...
        'sbatch_time',5,'sbatch_memory',16,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);        
    %get job id
    temp_job_id = erase(response.command_result{1},'Submitted batch job ');
    
%     %Fit motifs and cross-validate (parallellize by chunk)
%     [swarm_id, swarm_motifs] = FitMotifs_SpockSwarm(file_list_processed{cur_file},temp_job_id,s_conn,parameter_class);
%     job_id{cur_file} = [swarm_id{:}]; 
%     file_list_motifs{cur_file} = swarm_motifs;
end

%if you want to see what the postprocessed data looks like then run
%InspectPreprocessedData(PreprocessedDataFilepath,'postprocessed')

%% Cluster Motifs

%% Refit Clustered Motifs to the data

%% Refit the motifs built from all the animals to the data
job_id = [];
% Spock_RefitBasisMotifs_Swarm(file_list_processed,'Z:\Rodent Data\Wide Field Microscopy\ExampleData\',job_id,s_conn,parameter_class);
Spock_RefitBasisMotifs_Swarm(file_list_processed,'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL',job_id,s_conn,parameter_class);


%%

%Run clean up script and close out connection
ssh2_close(s_conn);
clear username password sconn




















