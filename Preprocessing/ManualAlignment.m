function [opts] = ManualAlignment(ref_img,opts)
ref_img = double(ref_img);

%show image
figure('name','Move rectangle and double-click to crop')
imagesc(ref_img); colormap gray

%Set to correct aspect ratio and shift to ~center of 1080p screen
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1)-200, pos(2)-200, size(ref_img,2), size(ref_img,1)]);
hold on

%cropping rectangle
rect = imrect(gca,[0 0 opts.crop_w opts.crop_h]);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(rect,fcn); 
setFixedAspectRatioMode(rect,1);
setResizable(rect,0);
crop_position = wait(rect);
close

% Crop the ref image
ref_img = imcrop(ref_img,crop_position);

%Get angle of brain 
angle = setMidlineAngle(ref_img);

%Roate the brain to vertical 
if angle <= 90 %if angled to right of yaxis
    alligned_img = imrotate(ref_img,-angle);
else %if angled left of yaxis
    alligned_img = imrotate(ref_img,180-angle);
end

%pad to keep cropping within bounds across lots of 
%different image sizes after rotating
alligned_img = padarray(alligned_img,[60,60]); 
bregma = setBregma(alligned_img);

%get coordinates for conservative cropping. 
%This make all recordings identically positioned with respect to bregma
x_crop_cord = (bregma(1)-opts.x_bregma_margin);
y_crop_cord = (bregma(2)-opts.y_bregma_margin); 

%crop the alligned image
cropped_alligned_img = double(imcrop(alligned_img,[x_crop_cord,y_crop_cord,...
    opts.crop_w,opts.crop_h]));

%crop can lead to inconsistent # of pixels. Make equal. 
cropped_alligned_img = cropped_alligned_img([1:opts.crop_h],[1:opts.crop_w]); 

%shift bregma to new coords in cropped image
bregma = [bregma(1)-x_crop_cord,bregma(2)-y_crop_cord];

%store variables in preprocessing_log structure
opts.bregma = bregma; 
opts.crop_cord = [x_crop_cord,y_crop_cord];
opts.alligned_img = alligned_img;
opts.cropped_alligned_img = cropped_alligned_img;
opts.angle = angle; 
opts.crop_position =crop_position; 


end





