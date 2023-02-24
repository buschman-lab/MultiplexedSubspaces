function fn=ParseBehavioralVideo(cur_rec,savedir)
%Camden MacDowell - timeless
%this is code for the 2021 experiments. 
%matches video to imaging frames and saves off different manual ROIs
if nargin <2; savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BehavioralData\'; end

[rec_name,~,~,~,~,FaceCam,~,tp] = LoadDataDirectories(cur_rec);

%manual inspection revealed that the behavioral video is not sampling
%specifically at 60 hz. So take the start and stop times and interpolate
%(see bottom for code). 
t = round(linspace(tp(1),tp(2),81000)); %81000 is the number of imaging frames

%adjust for the deconvolution
t = t(15+1:end-15);

%match to subspace downsample + 1 since derivative
t = t(1:2:end);

%+1 sample since derivative
t = [t,t(end)+median(diff(t))];

%load the video 
v = VideoReader(FaceCam);

%read the middle frame (so timing signal will already be on at this point. 
ref_img = read(v,t(1));

roi_names = {'face','whiskpad','nose','shoulder'};
% roi_names = {'whiskpad','nose','shoulder'};
%select the rois. 
roi = cellfun(@(x) SelectROI(ref_img,x), roi_names,'UniformOutput',0);
   
%show all rois
figure; imagesc(ref_img); 
cellfun(@(x) rectangle('Position',x.position,'EdgeColor',x.color,'FaceColor',[x.color 0.2],'LineWidth',2),roi,'UniformOutput',0);
axis off
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},[rec_name,'_ROI'],savedir,0); close all

%break frames into multiple 5 min blocks 
%5 min = ~6000 frames
f = 2000;
t_temp = t;
chunk = cell(1,floor(numel(t_temp)/f)+1);
for i = 1:floor(numel(t)/f)+1
    if i==floor(numel(t)/f)+1
        chunk{i} = t_temp; 
    else
        %add one timepoint to the end of each chunk since we are taking the derivative
        chunk{i} = t_temp(1:f+1);
        t_temp(1:f)=[];
    end
end

%save off data 
fn = [savedir,rec_name,'_parsedVideo.mat'];
save(fn,'t','chunk','FaceCam','roi_names','roi');




end %function 



% %quick snippet to inspect video and find start and stop time. 
% v = VideoReader(BodyCam);
% figure; 
% for i = 340384:340484
%    frame = read(v,i);
%    imagesc(frame);
%    title(sprintf('%d',i));
%    pause();
% end