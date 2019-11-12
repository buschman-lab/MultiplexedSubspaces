%add paths
addpath(genpath('C:\Users\Camden\Documents\Github\MotionMapper'));
addpath(genpath('C:\Users\Camden\Documents\Github\Widefield_Imaging_Analysis'));

%set filepaths
fn_path = 'C:\Users\Camden\Desktop\Videos\Mouse9031_10_17\';
fn_facecam = 'Cam_0_20191017-155226.avi';
fn_dlc = 'Cam_1_20191017-155226_Mouse431_10_17_2019DLC_resnet50_Headfixed_Behavior_BodyNov8shuffle1_120000.csv';
fn_savebase = 'Mouse9031_10_17';

%load behavioral analysis paramters
bp = behavioral_params; 

%parse the facecam to get timing signal and behavioral features
[data_mean, data] = ParseVideos([fn_path, fn_facecam],bp);

%save the roi image
save


%load dlc data


%filter, trim


%add the whiskertrace to the dlc data


%run wavelet transform and the tsne embedding


%watershed of phenograph and plot


%load motif H weightings (optionally scale by weightings?)


%find the onsets and take the 10 sec window around it


%for each motif, plot these trajectories in tnse space. 


%phenograph to cluster these trajectories per motif


%get the number of behavioral states for each motif 


