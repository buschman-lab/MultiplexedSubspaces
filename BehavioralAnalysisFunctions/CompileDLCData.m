function [features_downsampled] = CompileDLCData(T,ID)

switch ID
    case '431-10-17'
        fn_path = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse431_10_17_2019\';
        fn_facecam = 'Cam_0_20191017-155226.avi';
        fn_bodycam = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
        fn_dlc = 'Mouse431_10_17_2019DLC - XLSX.xlsx';
        fn_savebase = 'Mouse431_10_17_2019';
    case '432-10-18'
        fn_path = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse432_10_18_2019\';
        fn_facecam = 'Cam_0_20191018-151602.avi';
        fn_bodycam = 'Cam_1_20191018-151602_Mouse432_10_18_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
        fn_dlc = 'Mouse432_10_18_2019DLC-XLSX.xlsx';
        fn_savebase = 'Mouse432_10_18_2019';   
    case '432-10-17'
        fn_path = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\DeepLabCut_BehavioralState_Analysis\Mouse432_10_17_2019\';
        fn_facecam = 'Cam_0_20191017-171729.avi';
        fn_bodycam = 'Cam_1_20191017-171729_Mouse432_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000_labeled.mp4';
        fn_dlc = 'Mouse432_10_17_2019DLC-XLSX.xlsx';
        fn_savebase = 'Mouse432_10_17_2019';   
end

bp = behavioral_params;

%% Preprocessing Behavioral Videos 
%parse the facecam to get timing signal and behavioral features
if exist([fn_path fn_savebase 'roi_info.mat'],'file')
    roivals = load([fn_path fn_savebase 'roi_info.mat']);
    [facecam_mean, facecam_data, ~] = ParseVideos([fn_path, fn_facecam],bp,[],roivals.roi);
    onset = roivals.onset;
    offset = roivals.offset;
else
    [facecam_mean, facecam_data, roi] = ParseVideos([fn_path, fn_facecam],bp);
    %save the roi image
    axis off
    handles = get(groot, 'Children');
    saveCurFigs(handles,'-svg','Chosen ROIs',fn_path,1);
    close all; 

    %manual theshold for timing trace
    [onset, offset] = ManualThreshold(facecam_mean{1});

    %save off the roi information 
    save([fn_path fn_savebase 'roi_info.mat'],'roi','onset','offset');
end

%get the average motion energy (absolute derivative) of the pixels
temp = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), facecam_data(2:end),'UniformOutput',0);
face_motion_energy = cellfun(@(x) mean(abs(diff(x,2)),1)', temp,'UniformOutput',0);

%% Processing DLC Analyzed Data
%load dlc traced data
[~,~,raw_data] = xlsread([fn_path, fn_dlc]);

%get the speed of the 4 limbs
[limbs, id, ~] = parse_dlc(raw_data,{'frontrightpawcenter','frontleftpawcenter','backrightpawcenter','backleftpawcenter'},[],bp.dlc_epsilon);
limb_speed = cellfun(@(x) [0; mean(abs(diff(limbs(:,strcmp(id,x)),1)),2)],{'frontrightpawcenter','frontleftpawcenter','backrightpawcenter','backleftpawcenter'},'UniformOutput',0); 
limb_speed = [limb_speed{:}];
limb_speed = (mean(limb_speed,2));

%combined all features
face_motion_energy{2} = -1 * face_motion_energy{2}+max(face_motion_energy{2}(:)); %may need to flip the whisker energy if high whisking actually blurs the camera and make low energy  
features = cat(2,face_motion_energy{:},limb_speed);
labels = {'nose motion energy','whisker motion energy','limb speed'};
labels_abbrev = {'NME','WME','LS'};

%Trim to match start and stop of imaging 
features = features(onset:offset,:);

% for i = 1:size(features,2)
%     features(:,i) = convn(features(:,i),ones(130,1)/130,'same');
% end

% Downsample to match motif duration
features_downsampled = NaN(T,size(features,2));
for i = 1:size(features_downsampled,2)
    temp = features(:,i);
    features_downsampled(:,i) = interp1(1:numel(temp), temp, linspace(1,numel(temp),T),'linear'); 
end

%store the mapping from original to downsampled 
x_query_ds = linspace(1,size(features,1),T);

%log transform limb_speed;
features_downsampled(:,3) = log(features_downsampled(:,3));
features_downsampled = zscore(features_downsampled,1);



%%


















