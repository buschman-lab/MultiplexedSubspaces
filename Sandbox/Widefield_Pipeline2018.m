

%select folders to process and grab the first file from each rec. 
[file_list,folder_list] = GrabFiles('Pos0.ome.tif');

%configure preprocessing options
opts = ConfigurePreProcessing();
 
%Grab reference images for each. Preload so no delay between loop.
ref_imgs = cellfun(@(x) GetReferenceImage(x,opts.fixed_image),...
    file_list, 'UniformOutput',0);

for cur_fold = 1:numel(folder_list)
    fprintf('\n Working on recording %d of %d',cur_fold, numel(folder_list));
    clear prepro_log
    %manual allignment 
    prepro_log = ManualAlignment(ref_imgs{cur_fold},opts);

    %mask vasculature and manual cleanup (optional)
    prepro_log = MaskVasculature(...
        prepro_log.cropped_alligned_img,prepro_log);
    
    %save off the options to each folder
    save([folder_list{cur_fold} filesep 'prepro_log'],'prepro_log')
end

%% Write a bash sript to send the recordings in for processing

%Preprocesses the recordigns

%load the opts file
%apply the cropping and allignment to every frame;
%downsample
%split wavelength
%get calculate dff 




%% Section 3: Preprocess Raw Image Stacks (Downsample, Crop, Allign, Make DFF, detrend)
AllOpts(1:length(file_list)) = struct('Opts',struct());
for cur_file = 1:length(file_list) %File processing loop
    %Load options file either from the designated save dir (if exists) or
    %from the image stack directory (if no designated save dir); 
    [path, base] = fileparts(file_list{cur_file});
    %remove the rec # denominator to load the same opts file for multiple
    %recs on same day
    base = erase(base,'_MMStack_Pos0.ome');
    base = base(1:end-2);
    if isempty(SaveDir)
        opts_fn = [path filesep sprintf('%s_Processing',base),filesep, 'opts.mat'];
    else
        opts_fn = [SaveDir filesep sprintf('%s_Processing',base), filesep, 'opts.mat'];
    end
    if ~exist(opts_fn,'file')        
        error('No options file for file %s',file_list{cur_file})
    else
        load(opts_fn)
    end
    AllOpts(cur_file).Opts = opts;
end

if ParProcess ==1 
    parfor cur_file = 1:length(file_list)
        tic
        %Get the specific name of each file (for days with multiple recs) so
        %they don't overwrite when saved as dffs
        [path, base] = fileparts(file_list{cur_file});
        AllOpts(cur_file).Opts.DFF_fn = [AllOpts(cur_file).Opts.SaveDir filesep sprintf('%sDFF.mat',base)];
        AllOpts(cur_file).Opts.Stack_fn = [AllOpts(cur_file).Opts.SaveDir filesep sprintf('%s_PreProStack.mat',base)];
        %If parralel processing, showing transformed info will stop program
        if AllOpts(cur_file).Opts.showTransformed
           AllOpts(cur_file).Opts.showTransformed = 0; 
        end
        Widefield_PreProcessRS(file_list{cur_file},AllOpts(cur_file).Opts);
        toc
        fprintf('Done processing %s... it took %d minutes',base,(toc./60));
    end %file processing loop
else
    for cur_file = 1:length(file_list)
        tic
        %Get the specific name of each file (for days with multiple recs) so
        %they don't overright
        [path, base] = fileparts(file_list{cur_file});
        AllOpts(cur_file).Opts.DFF_fn = [AllOpts(cur_file).Opts.SaveDir filesep sprintf('%sDFF.mat',base)];
        %If parralel processing, showing transformed info will stop program
        if AllOpts(cur_file).Opts.showTransformed
           AllOpts(cur_file).Opts.showTransformed = 0; 
        end
        Widefield_PreProcessRS(file_list{cur_file},AllOpts(cur_file).Opts);
        toc
        fprintf('Done processing %s... it took %d minutes',base,(toc./60));
    end %file processing loop
end
%%

















