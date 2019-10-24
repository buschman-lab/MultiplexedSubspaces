%Initializes everything for the StimulusEvoked_Pipeline2019_Spock pipeline
%Define options
SaveDir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Processed_SensoryEvoked'; %set the save directory

pp = genpath('Z:\Rodent Data\Wide Field Microscopy\AnalysisCode_Repository\');
addpath(pp);

% load([SaveDir filesep 'CompleteFile_List.mat']);

cd('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\SensoryEvokedMapping');

%Select Imaging Data To Process
tempd = dir(pwd);
isub = [tempd(:).isdir]; %returns logical vector
nameFolds = {tempd(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = []; %remove silly extra folders

%Sensory evoked is in sub-sub-directories so move into each one
file_list = [];
for cur_dir = 1:numel(nameFolds)
   fprintf('\n Working through folder %d of %d',cur_dir,numel(nameFolds));
   tempd = dir(nameFolds{cur_dir});
   isub = [tempd(:).isdir]; %returns logical vector
   targetdir = {tempd(isub).name}';
   targetdir(ismember(targetdir,{'.','..'})) = []; %remove silly extra folders
    
   for cur_targ = 1:numel(targetdir)
      [temp_file_list] = SelectImagingData([nameFolds{cur_dir} filesep targetdir{cur_targ}]);
      file_list = [file_list temp_file_list];
   end
end
file_list = unique(file_list); 
%Just get file names to rename them for spock processing
bucket_fn = [];
for i = 1:numel(file_list)
    file_list{i}(1:3) = [];
    file_list{i} = strrep(file_list{i},'\','/');
    bucket_fn{i} = ['/jukebox/buschman/', file_list{i}];  
end
save([SaveDir filesep 'CompleteFile_List.mat'],'file_list','bucket_fn');

%Initialize options and define ROI on the data
for cur_file = 1:length(file_list) %File Processing Loop
    Initialize_Sensory(file_list{cur_file},SaveDir); 
    fprintf('Initialized Preprocessing info for file %d, out of %d...\n',cur_file,length(file_list))
end %End File Processing Loop



%% code to double check the brains
% %move to folder of interest
% figure;
% cd(SaveDir);
% tempd = dir(pwd);
% isub = [tempd(:).isdir]; %returns logical vector
% nameFolds = {tempd(isub).name}';
% nameFolds(ismember(nameFolds,{'.','..'})) = []; %remove silly extra folders
% COUNT = 1;
% for i = 1:numel(nameFolds)
%    cd(SaveDir)
%    cd(nameFolds{i});
%    load opts
%    subplot(10,8,COUNT);COUNT = COUNT+1;
%    imagesc(opts.CroppedAllignedImg);
%    title(sprintf('%s',opts.RecordingName));
%    axis square
% end

















