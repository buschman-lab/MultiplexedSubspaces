function initializePrePro(in_fns,SaveDir)
%Camden MacDowell 2018
%% Build options structure
close all
if nargin <2
    opts.SaveDir ='';
else
    opts.SaveDir = SaveDir;
end
%INITIALIZATION OPTIONS
opts.Verbose = 1; %whether to be chatty about what we are doing - default to 1 
opts.OverwriteInitialization = 1; %overwrite previous initialiation?
opts.CropHeight = 540; %height of croped image in pixels. 
opts.CropWidth = 540; %width of cropped image in pixels
opts.FixedImage =1; %What to use for the ref image? 1 = middle image (default), 0 = meanprojection (slow). 
opts.MaskVasculature = 1; %Mask vasculature
opts.ManualMask = 1; %Add additional manual masking to image (e.g. if imperfection in skull); 
opts.VascSTD = 2.5; %number of standard devidations from mean of ref_img to consider vasculature (2.5 works well)
opts.CloseDiskSize = 2; %pixel size of disk used to remove salt/pepper noise on vascmask. 
opts.MaskBrainOutline = 1; %Create mask for outside of brain (manually made). 
opts.MaskBrainOutlineDir = 'C:\Users\macdo\OneDrive\Buschman Lab\AnalysisCode_Repository\PreProcessing Code\brainoutline.mat';
% opts.MaskBrainOutlineDir = 'C:\Users\Camden\Documents\Work\OneDrive\Buschman Lab\AnalysisCode_Repository\PreProcessing Code\brainoutline.mat';
opts.XBregmaMargin = 280;  %For outlininign brain: Number of pixels to anterior to bregma (default = 280)
opts.YBregmaMargin = 230; %for outlining brain: number of pixels to keep lateral to bregma (default - 230)
%Save directory info 
if isempty(opts.SaveDir) %if no specified dir, then save to imagestack location
    [cur_pn, opts.SaveFileBase, ~] = fileparts(in_fns);
    opts.SaveDir = [cur_pn filesep sprintf('%s_Processing',opts.SaveFileBase)];
else
    [~, opts.SaveFileBase, ~] = fileparts(in_fns);
    opts.SaveDir = [opts.SaveDir filesep sprintf('%s_Processing',opts.SaveFileBase)];
end

%PREPROCESSING OPTS
opts.makeDFF = 1; %Should we make a dff of the recordings (this will make a dff of whole rec and not parse by motion/nonmotion epochs
opts.dsSpatial = 1; %should we spatially downsample?
opts.dsFactor = 4; %Amount to ds(if ends @128x128 ~=80um^2/px)
opts.motion = 0; %Do motion correction (with image registration)
opts.UsePhaseRegistration = 1; %if true, then uses imregcorr, otherwise uses imregister
opts.TransformationType = 'similarity'; %can be 'similarity', 'rigid', or 'translation'
opts.gpu = 0; %Perform movement correction on the gpu using imregdemons (displacement field)
opts.showTransformed = 0; 
opts.crop = 1;%1; %Crop during PreProcessing 
opts.PreProVerbose = 1; %whether to be chatty about what we are doing
opts.method = 'movingavg'; %'mean','median','mode', 'movingavg' for dff: moving avg is set to calculate a 10sec sliding window mean baseline 
opts.FrameStop = 0; %stop calculating the dff at a particular frame
opts.Detrend = 1; %linear detrending of the dff (removes slow flucts in signal due to LED dynamics, etc.)
opts.CalcBaselineForSpecificIndx = 0; %Only calc baseline f along a praticular frame index
opts.BaselineIndx = [1,30]; %the above indx
opts.ApplyMask = 1;
opts.dsIgnoreNan = 1; % 1= Ignore vasc mask Nan when downsampling (e.g. nansum vs sum)
opts.Fs = 13;
opts.OnlyInitializeOptsForFirstRec = 1;

if opts.MaskBrainOutline %load brainoutline mask
    %load a premade 'outline' logical mask of brain. 1 = brain 0 = other
    load(opts.MaskBrainOutlineDir);   
end

%Only initialize for the first rec of the day;
if opts.OnlyInitializeOptsForFirstRec
    IsFirstRec = strfind(opts.SaveFileBase,'_1_MMStack_Pos0.ome'); 
    if isempty(IsFirstRec) %exit initialization if this isn't the first rec of the day
        warning('Already have opts file (from epoch 1)... exiting initialization'); 
    return
    end
end

%Remove the file part deliniator from the opts save dir so all recs from a day are grouped
%together
opts.SaveDir = erase(opts.SaveDir,'_1_MMStack_Pos0.ome');
opts.RecordingName = erase(opts.SaveFileBase,'_1_MMStack_Pos0.ome');

opts.Opts_fn = [opts.SaveDir filesep sprintf('opts.mat')];
%check to see if options already exist
if exist(opts.Opts_fn)
    if opts.OverwriteInitialization
    else
        warning('Using Previous Initialization Data for data %s',opts.SaveFileBase); 
        return
    end
end

if ~exist(opts.SaveDir, 'dir')
    mkdir(opts.SaveDir);
end


%% Count the number of image in the input file, grab the ref image and process
%Surpress warning regarding tag structure of tiff. 
warning ('off','all');
in_tiff = Tiff(in_fns, 'r');

%Count how many images in each file
in_tiff.setDirectory(1);
while ~in_tiff.lastDirectory
    in_tiff.nextDirectory;
end
img_count = in_tiff.currentDirectory;

%Grab the middle image
if opts.FixedImage == 1
    img_idx = floor(sum(img_count)/2);

    fprintf('Using middle image as ref image (%d of %d)\n', img_idx, img_count);
    cur_in_file = find(img_idx <= cumsum(img_count), 1, 'first');
    if cur_in_file > 1
        in_tiff(cur_in_file).setDirectory(img_idx - sum(img_count(1:(cur_in_file-1))));
    else
        in_tiff(cur_in_file).setDirectory(img_idx);
    end
    ref_img = in_tiff(cur_in_file).read;
elseif opts.FixedImage ==0
    fprintf('Calculating average image of first % to use as fixed image...\n');
    fixed_img = zeros(in_tiff(1).getTag('ImageLength'),in_tiff(1).getTag('ImageWidth'));
    stack = uint16(zeros(size(fixed_img,1),size(fixed_img,2),img_count));
    for cur_img_ind = 1:img_count
        %Find the current input image file
        cur_in_file = find(cur_img_ind <= cumsum(img_count), 1, 'first');
        if cur_in_file > 1
            in_tiff(cur_in_file).setDirectory(cur_img_ind - sum(img_count(1:(cur_in_file-1))));
        else
            in_tiff(cur_in_file).setDirectory(cur_img_ind);
        end
        %Read the image
        cur_img = double(in_tiff(cur_in_file).read);
        %Save to overall variable
        try
            fixed_img = fixed_img(:,:,1) + cur_img(:,:,1);
        catch
            error('Dimensions issue combining Fixed_img and Cur_img. Error occurs on %g',cur_img_ind);
        end
        if mod(cur_img_ind, round(0.1*sum(img_count))) == 0
            fprintf('\t%2.0f%% done...\n', (cur_img_ind./sum(img_count)*100));
        end
        stack(:,:,cur_in_file) = (in_tiff(cur_in_file).read); 
    end %image loop
    if any(fixed_img >= (2^64 - 1)), error('Too many images to average at once...'); end
    %Calculate and write average image
    ref_img = fixed_img./sum(img_count);
end
%turn on warnings again
warning ('on','all');

%Quick check to make sure ref_img is not smaller than 512x512 
if size(ref_img)<=512
    warning('%s is too below 512x512 which is incompatible with this processing pipeline',opts.SaveFileBase);
end
%% Set cropping, rotation, and bregma information;
ref_img = double(ref_img);
figure('name','Move rectangle and double-click to crop')
imagesc(ref_img); colormap gray
%Set to correct aspect ratio and shift to ~center of 1080p screen
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1)-200, pos(2)-200, size(ref_img,2), size(ref_img,1)]);
hold on
title(sprintf('Crop %s...',opts.SaveFileBase));
rect = imrect(gca,[0 0 opts.CropHeight opts.CropWidth]);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(rect,fcn); 
setFixedAspectRatioMode(rect,1);
setResizable(rect,0);
crop_position = wait(rect);
close
% Crop the ref image
ref_img = imcrop(ref_img,crop_position);
%Store Ref_img and crop info
opts.ref_img = double(ref_img);
opts.crop_position= crop_position;
%Get angle of brain 
angle = setMidlineAngle(ref_img);
%Keep angle info
opts.angle = angle; 
%Roate the brain to vertical 
if angle <= 90 %if angled to right of yaxis
    AllignedImg = imrotate(ref_img,-angle);
else %if angled left of yaxis
    AllignedImg = imrotate(ref_img,180-angle);
end

%pad to keep cropping within bounds across lots of different iamge sizes after rotating
AllignedImg = padarray(AllignedImg,[60,60]); 
Bregma = setBregma(AllignedImg);
%get coordinates for conservative cropping (e.g. leaving space on side of brain hence 540x540).
opts.XCropCord = (Bregma(1)-opts.XBregmaMargin); %280
opts.YCropCord = (Bregma(2)-opts.YBregmaMargin); %230
CroppedAllignedImg = double(imcrop(AllignedImg,[opts.XCropCord,opts.YCropCord,opts.CropWidth,opts.CropHeight]));
CroppedAllignedImg = CroppedAllignedImg([1:opts.CropWidth],[1:opts.CropHeight]); 
%shift bregma to new coords in cropped image
Bregma = [Bregma(1)-opts.XCropCord,Bregma(2)-opts.YCropCord];
opts.Bregma = Bregma; %Store Bregma
opts.CroppedAllignedImg = CroppedAllignedImg; %Store Alligned Image

%% Identify Vasculature
if opts.MaskVasculature
    %Indentify vasculature by subtracting a large median filter  (e.g. smoothed image
    %form the base image. This excentuates dark areas that are surrounded by
    %light areas (e.g. vasculature). 
    smoothed = medfilt2(CroppedAllignedImg,[125 125]); %Use large neighboorhood e.g. 125; 
    Mask = CroppedAllignedImg-smoothed; 

    %threshold to remove pxs lower that 'x' standard deviations below mean
    %of VascMask (e.g. the vasculature and other dark blemishes). 
    mVal = nanmean(nanmean(Mask));
    sVal = nanstd(nanstd(Mask)); 
    Mask = Mask>=(mVal-(opts.VascSTD*sVal));

    %close the image to clean up noise
    se = strel('disk',opts.CloseDiskSize);
    Mask = imclose(Mask,se); 

    %Add a premade mask to outline the brain.
    if opts.MaskBrainOutline 
        Mask = Mask+brainoutline;
        Mask = Mask==2;
    end

    if opts.Verbose || opts.ManualMask
        imshowpair(CroppedAllignedImg,CroppedAllignedImg.*Mask);
        hold on; 
        scatter(Bregma(1),Bregma(2),'*');
        title('Cropped Img with VascMask overlay');
    end

    %Manually Nan regions and unmask others 
    if opts.ManualMask
        dlg = questdlg('Would you like to manually mask regions?', ...
        'MASK ROIS','Yes','No','No');
        while 1
            switch dlg
                case 'Yes' %Mask additional region
                    close all                
                    figure('name','Draw ROI around additional regions to NaN');
                    imagesc(CroppedAllignedImg.*Mask)
                    ManualMask = impoly;
                    wait(ManualMask);               
                    ManualMask = ManualMask.createMask();
                    Mask(ManualMask)=0; close; 
                    figure('name','Current Masked Image')
                    imshowpair(CroppedAllignedImg,CroppedAllignedImg.*Mask); 
               case 'No'; break %end dlg
            end
            choice = questdlg(sprintf('Would you like to mask additional regions'),...
                'MASK ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
        dlg = questdlg('Would you like to manually UNmask regions?', ...
        'MASK ROIS','Yes','No','No');
        while 1 %Manuall UNmask regions 
            switch dlg
                case 'Yes' %Mask additional region
                    close all
                    figure('name','Draw ROI around additional regions to NaN');
                    imagesc(CroppedAllignedImg.*Mask)
                    ManualMask = impoly;
                    wait(ManualMask);               
                    ManualMask = ManualMask.createMask(); %his is 1 in the desired ROI          
                    Mask(ManualMask)=1; close
                    figure('name','Current Masked Image')
                    imshowpair(CroppedAllignedImg,CroppedAllignedImg.*Mask);
               case 'No'; break
            end
            choice = questdlg(sprintf('Would you like to unmask additional regions'),...
                'MASK ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
    end

    if opts.Verbose || opts.ManualMask
        close all
        imshowpair(CroppedAllignedImg,CroppedAllignedImg.*Mask);
        hold on; 
        scatter(Bregma(1),Bregma(2),'*');
        title('FINAL cropped Img with VascMask overlay');
            choice = questdlg(sprintf('Would you like to flag this brain for imperfections?'),...
            'Flag brain',...
            'Yes','No','No');
       switch choice
           case 'Yes'
               opts.flaggedBrain = 1;
           case 'No' 
               opts.flaggedBrain = 0; 
      end
    end
    %Store mask
    opts.Mask = Mask; 
    %Store the complete masked image for reference
    opts.MaskedRefImg = CroppedAllignedImg.*Mask;
elseif opts.MaskVasculature==0 && opts.MaskBrainOutline == 1 %only apply exterior mask
    Mask = brainoutline;
    if opts.Verbose
        close all
        imshowpair(CroppedAllignedImg,CroppedAllignedImg.*Mask);
        hold on; 
        scatter(Bregma(1),Bregma(2),'*');
        title('FINAL cropped Img with Masked Brain Overlay');
    end
    %Store mask
    opts.Mask = Mask; 
    %Store the complete masked image for reference
    opts.MaskedRefImg = CroppedAllignedImg.*Mask;
else %apply no mask
    Mask = [];
    opts.MaskedRefImg = [];
end
    

%Save off options
save(opts.Opts_fn,'opts');
fprintf('\tDone initializing imagestack %s...\n',opts.SaveFileBase)

end








