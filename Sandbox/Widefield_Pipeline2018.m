

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
    
    %close figure
    close
    
    %save off the options to each folder
    save([folder_list{cur_fold} filesep 'prepro_log'],'prepro_log')
end

%% Call spock function to Loop through the files and process each recording
%get a all the recordings in each folder
[file_list,folder_list] = GrabFiles('.tif',0,folder_list{cur_fold});

%load opts function;

stack = PreProcess(in_fn,opts);
if numel(unique(opts.wavelength_pattern))>1 %if multiple wavelengths used
   dff = HemodynamicCorrection(stack, opts); 
else
   dff = makeDFF(stack, opts);
end
    
%Save off the dff with pre-pended file name





