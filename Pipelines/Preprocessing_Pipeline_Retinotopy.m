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
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis\'));

%configure preprocessing options
opts = ConfigurePreProcessing('crop_w',540,'vasc_std',1,'save_uncorrected',1,'fixed_image','first','method','baseline');

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_retinotopy';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

%set up save directories
save_dir_processed = 'Z:\Projects\Cortical Dynamics\Retinotopy\Processed\'; %target savedirector
if ~exist(save_dir_processed,'dir')
    mkdir(save_dir_processed);
end


%% Manual Portion 
rec_path = 'Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\Retinotopy\Mouse360_10_14_2021';
warning('can only be used on one recording (i.e. one set of folders) as a time');
%select folders to process and grab the first file from each rec.
%EXAMPLE DATA: Select 'folders' and then select 'Z:\Rodent Data\WideField Microscopy\ExampleData\Mouse431_10_17_2019\431-10-17-2019_1'
[file_list_first_stack,folder_list_raw] = GrabFiles('Pos0.ome.tif',1,{rec_path});

%Grab reference images for each. Preload so no delay between loop.
ref_img = GetReferenceImage(file_list_first_stack{1},opts.fixed_image);

%loop through reference images, register, and apply manual changes
%manual allignment 
prepro_log = ManualAlignment(ref_img,opts);

%mask vasculature and manual cleanup (optional)
prepro_log = MaskVasculature(...
    prepro_log.cropped_alligned_img,prepro_log);

close
%no transformation
prepro_log.tform = []; 
prepro_log.output_size = [];

%grab the trial timing and parse and save off for spock. these timings are
%the same across trials (assumes no dropped frames, which is a safe
%assumption in my experience)
[LogFn,log_folder] = GrabFiles('acquisitionlog.m',0,{rec_path});
prepro_log.frame_index_blue = parseRetinoTiming(LogFn{1},'removeHemo',1);
prepro_log.frame_index = parseRetinoTiming(LogFn{1},'removeHemo',0);

%stimulus type info
[LogFn,~] = GrabFiles('stimInfo.mat',0,log_folder);
temp = load(LogFn{1});
prepro_log.stimInfo = temp.stim_type;

%save off the options to overarching folder
save([log_folder{1} filesep 'prepro_log'],'prepro_log')

%% Preprocessing. Results in a single hemo corrected, masked recording for each day in the 'preprocessed' folder

%run as a loop of spock jobs
[opts_path,~] = GrabFiles('prepro_log.m',0,log_folder); 
folder_trials = dir(log_folder{1});
folder_trials = folder_trials([folder_trials(:).isdir]);
folder_trials = folder_trials(~ismember({folder_trials(:).name},{'.','..'}));
folder_trials = arrayfun(@(n) [log_folder{1},filesep,folder_trials(n).name],1:size(folder_trials,1),'UniformOutput',0);

for cur_trial = 1:numel(folder_trials) %trials in separate folders    
   [file_list_raw,~] = GrabFiles('.tif',0,folder_trials(cur_trial));
   input_val = {ConvertToBucketPath(file_list_raw{1}), ConvertToBucketPath(opts_path{1})};
   script_name = WriteBashScript(parameter_class,sprintf('%d_%d',1,cur_trial),'Spock_Retinotopy_Pipeline',input_val,{"'%s'","'%s'",'%d'},...
       'sbatch_time',5,'sbatch_memory',12);  %

   %Run job
   response = ssh2_command(s_conn,...
       ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
       sprintf('sbatch %s',script_name)]); 

   %get job id
   job_id{cur_trial} = erase(response.command_result{1},'Submitted batch job ');
   if cur_trial ~=numel(file_list_raw)
       job_id{cur_trial} = [job_id{cur_trial} ','];
   end   
end

%% (need to write spock for this)
%load all the files
fn = arrayfun(@(x) GrabFiles('stack.mat',0,x),folder_trials,'UniformOutput',0);

%load the data 
index = cellfun(@(x) load(x,'index'),[fn{:}],'UniformOutput',1);
index = [index(:).index];

%stimulus ids
opts = load(fn{1}{1},'opts');
stimInfo = opts.opts.stimInfo(index); %sorted by the index

%load data
data = cellfun(@(x) load(x,'dff'),[fn{:}],'UniformOutput',0);
data = cellfun(@(x) x.dff,data,'UniformOutput',0);

%truncate to same length
minlength = min(cellfun(@(x) size(x,3),data,'UniformOutput',1));

%data 
data = cellfun(@(x) x(:,:,1:minlength),data,'UniformOutput',0);

%combine by stim type
data_avg = arrayfun(@(x) nanmean(cat(4,data{stimInfo==x}),4),unique(stimInfo),'UniformOutput',0);

%% extract retinotopy positions for each pixel
Fs = 15;
T=1/Fs;
horz_degree = 60/5; %degree per sec
vert_degree = 40/5; %degree per sec

phase_img = zeros(68,68,4);
for j = 1:4
    [temp,nanpxs] = conditionDffMat(data_avg{j});
    temp=temp';
    temp_phase = [];
    for i = 1:size(temp,1)
        x = temp(i,:);
        x=x-mean(x);
%         x = [zeros(size(x)),x,zeros(size(x))];
        t_domain = (0:length(x)-1)/Fs;
        xw = 2*fft(x)/length(x);
        xw = xw(1:floor(length(xw)/2));
        xw_mag = abs(xw);
        xw_phase = angle(xw)*180/pi;
        f_domain = 0:(Fs/2)/(length(xw)-1):Fs/2;
        if j<=2
            temp_phase(i,:) = nanmean(xw_phase(2));
        else
            temp_phase(i,:) = nanmean(xw_phase(2));
        end
    end
    phase_img(:,:,j) = conditionDffMat(temp_phase',nanpxs);
end

%subtract opposite sides
phase_img = cat(3,phase_img(:,:,1)-phase_img(:,:,2),phase_img(:,:,3)-phase_img(:,:,4));
% phase_img = cat(3,phase_img(:,:,1),phase_img(:,:,3));

% Compute the visual field sign map
[dhdx, dhdy] = gradient(phase_img(:,:,1));
[dvdx, dvdy] = gradient(phase_img(:,:,2));

graddir_hor = atan2(dhdy,dhdx);
graddir_vert = atan2(dvdy,dvdx);

vdiff = exp(1i*graddir_hor) .* exp(-1i*graddir_vert); %Should be vert-hor, but the gradient in Matlab for y is opposite.
VFS = sin(angle(vdiff)); %Visual field sign map
id = find(isnan(VFS));
VFS(id) = 0;

hh = fspecial('gaussian',size(VFS),3);
hh = hh/sum(hh(:));
VFS = ifft2(fft2(VFS).*abs(fft2(hh)));  %Important to smooth before thresholding below
%%
[~,name] = fileparts(fn{1});
name = erase(name,'_1_stack');
% load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\MesoscaleDynamics_2020_repo_depreciated_with_redundancies\Code_From_Local_Machine\FigureMask.mat')
temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis\Preprocessing\brainoutline.mat');
mask = imresize(temp.brainoutline,[68 68]);
VFS = VFS(:);
VFS_three=VFS;
% mask(:,1:35)=0;
VFS(mask(:)==0)=NaN;
VFS = reshape(VFS,[68 68]);
% close; imagesc(VFS([40:65],[35:65]))
close; imagesc(-1*VFS,[-0.15 .15]);  
axis square
axis off
colormap magma
%overlay distance to bregma
b = floor(opts.opts.bregma/opts.opts.spatial_bin_factor);
hold on; 
plot(b(1)-0.5,b(2),'color','red','markersize',40,'marker','.')
%each pixel is ~135uM so overlay a scale line to where you are targeting
line([5,5],[5,5+7.4],'color','red','linewidth',3)
%%
%for V1
saveas(gcf,[log_folder{1}, filesep,name,'visualfieldsign.png']); close all; 
saveas(gcf,[log_folder{1}, filesep,name,'visualfieldsign.fig']); close all; 
% %saveoffthereference image as well
temp = opts.opts.cropped_alligned_img;
imagesc(temp,[prctile(temp(:),2.5),prctile(temp(:),97.5)]); colormap gray
axis off; axis square
saveas(gcf,[log_folder{1}, filesep,name,'ReferenceImage.png']); close all; 
%%
%Run clean up script and close out connection
ssh2_close(s_conn);
clear username password sconn














