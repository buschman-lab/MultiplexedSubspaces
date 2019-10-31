function [data_mean, data] = ParseVideos(in_fn)

bp = behavioral_params; 
%load video 
vid = VideoReader(in_fn);

%get the camera id, face vs body = 1 or 2
camera_idx = str2num(cell2mat(regexp(in_fn,'(?<=Cam_)\d+(?=_)','match')))+1; %+1 since 0 indexed 

%read the first frame
ref_img = read(vid,1);

%select the rois. 
roi = cellfun(@(x) SelectROI(ref_img,x), bp.roi_names{camera_idx},'UniformOutput',0);

%show all rois
figure; imagesc(ref_img); 
cellfun(@(x) rectangle('Position',x.position,'EdgeColor',x.color,'FaceColor',[x.color 0.2],'LineWidth',2),roi,'UniformOutput',0);

%Parse: framewise for possible memory issues
img_count = vid.Duration*vid.FrameRate; 

%Preallocate data structure
data = cellfun(@(x) single(NaN(x.position(4)+1,x.position(3)+1,img_count)),roi,'UniformOutput',0); %add one pixel to the dimensions for cropping
fprintf('\nParsing Video...')
for cur_img_ind = 1:img_count
    if mod(cur_img_ind,round(0.10*img_count)) ==0
        fprintf('\t%g%% Complete\n', round(cur_img_ind./img_count*100,2));
    end
    img = read(vid,cur_img_ind);
    for cur_roi = 1:numel(roi)
        if cur_roi==1 %use just the green channel for timing signal
            data{cur_roi}(:,:,cur_img_ind) = single(imcrop(img(:,:,2),roi{cur_roi}.position)); %use the second wavelength,it picks up the blue light best
        else
            data{cur_roi}(:,:,cur_img_ind) = single(rgb2gray(imcrop(img,roi{cur_roi}.position))); %all wavelengths weighted equally. 
        end
    end  
end

%Get mean trace per roi
data_mean = cellfun(@(x) nanmean(x,1), data,'UniformOutput',0);


end
%vectorize
%data = cellfun(@(x) reshape(x,[size(x,1)*size(x,2)],size(x,3)),data,'UniformOutput',0);

%%PENDING
%First x PCs for each roi 

    






















