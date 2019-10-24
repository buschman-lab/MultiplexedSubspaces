function StimulusEvoked_Pipeline2019_Spock(jobID)

%This does all preprocessing on spock
%requires that you first initialize (select bregma/etc.) for each file

bucket = '/jukebox/buschman/';

savedir = [bucket '/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/Processed_SensoryEvoked']; %set the save directory

if exist(savedir)==0
    mkdir(savedir);
end
fprintf('Working on Job ID %d',jobID);

%Add path of the entire analysis code directory 
pp = genpath([bucket 'Rodent Data/Wide Field Microscopy/AnalysisCode_Repository']);
addpath(pp);

%Load the file list of 
file_list = load([savedir filesep 'CompleteFile_List.mat']); 
file_list = file_list.bucket_fn;

%Divide the file_list up into 1000 jobs
x = (1:1:numel(file_list));
cur_grp = [];
COUNT = 1; 
for i = 1:1000
    temp = [];
    for ii = 1:ceil(numel(file_list)/1000)
        if COUNT>numel(file_list)
            break
        end
        temp(ii) = x(COUNT);
        COUNT = COUNT+1;
    end
    cur_grp{i} = temp;
end
cur_grp(cellfun('isempty',cur_grp)) = [];

if jobID>numel(cur_grp)
    fprintf('Job Array ID Exceeds cur_grp bounds')
    return
end

%% Preprocess the data
for cur_file = cur_grp{jobID}
    tic
    %Load options file either from the designated save dir (if exists) or
    %from the image stack directory (if no designated save dir); 
    [path, base] = fileparts(file_list{cur_file});
    DFF_fn = sprintf('%sDFF.mat',base);
    Stack_fn = sprintf('%s_PreProStack.mat',base);
    %remove the rec # denominator to load the same opts file for multiple
    %recs on same day
    base = erase(base,'_MMStack_Pos0.ome');
    base = base(1:end-2);
    %this is sloppy, but remove the '_' from the last line of the name if
    %it's there
    if strcmp(base(end),'_')
        base = base(1:end-1);
    end
    fprintf('This is the base\n %s \n',base);
    fprintf('This is the path\n %s \n',path);

    Opts_fn = ([savedir filesep sprintf('%s_Processing',base) filesep 'opts.mat']);
    savedir_dff = [savedir filesep sprintf('%s_Processing',base)];

    if ~exist(Opts_fn,'file')        
        error('No options file for file %s',file_list{cur_file})
    else
        opts = load(Opts_fn);
        opts = opts.opts;
    end
    %Get the specific name of each file (for days with multiple recs) so
    %they don't overright
    opts.showTransformed = 0;
    
    opts.BaselineIndx = [1 33];
    opts.Detrend =0;
    opts.DFF_fn = [savedir_dff filesep DFF_fn];
    opts.Stack_fn = [savedir_dff filesep Stack_fn];
    %If parralel processing, showing transformed info will stop program  
    Widefield_PreProcess_Sensory(file_list{cur_file},opts);
    toc
    fprintf('Done processing %s... it took %d minutes',base,(toc./60));
end
end
























