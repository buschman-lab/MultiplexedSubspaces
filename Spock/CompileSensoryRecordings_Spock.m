function CompileSensoryRecordings_Spock(file_path)

if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis/'));
end

[folder, fn] = fileparts(file_path); 
MouseNum = fn(6:9);

fn = erase(fn,'_MMStack_Pos0.omeDFF');
%trial number is between the final _ and the end
tempidx = regexp(fn,'_');
tempidx = tempidx(end); %only grab the final _ in case there are multiple
TrialNum = str2num(fn(tempidx+1:end));

%Find the stimulus type of the current trial
cd(folder);
opts_fn = dir('*OptsFile.mat');
opts_info = load(opts_fn.name);
StimType = opts_info.stimOrder(TrialNum);

%timing info
timing_fn = dir('*TimingFile.mat');
timing_info = load(timing_fn.name);

%Get the StimStartFrame 
StimStartF = floor(timing_info.timingdata.timetostim(TrialNum)*opts_info.opts.framerate);

data = load(file_path);
data = data.dff;

mask = load('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/Spock_Code_Repository/Functions/FigureMask_Vis.mat');
mask = mask.mask;

% Take window after sensory delivery
data = data(:,:,StimStartF-13:end);

%Downsample
data = SpatialBin(data,2);

%Apply mask
for i = 1:size(data,3)
    temp = data(:,:,i);
    temp(mask==0)=0;
    data(:,:,i) = temp;
end

%deconvolve and normalize
[data, nanpxs] = ProcessAndSplitData(data,[],'general_params_sensoryResponses');    

%save off
cd('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/VPA_Mesomapping/SensoryEvoked');
data = data';
save([fn,'_',num2str(StimType),'.mat'],'data','nanpxs','StimType','TrialNum','MouseNum','file_path','opts_info','timing_info');
fprintf('\nDONE');
end