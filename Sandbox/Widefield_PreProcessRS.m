function Widefield_PreProcessRS(in_fns,opts)


%% Loop through images, read each one, and perform the desired compuations and save
%Preallocate stack
if opts.ApplyMask == 1
    stack = zeros(size(opts.Mask,1),size(opts.Mask,2));
else
    if opts.crop
        stack = zeros(floor(opts.crop_position(3)),floor(opts.crop_position(4)));
    else
        stack = zeros(floor(size(ref_img,1)),floor(size(ref_img,2)));
    end
end
%Preallocate along the 3rd dimension
if opts.dsSpatial
    stack = SpatialBin(stack,opts.dsFactor);
else
    stack = repmat(stack, 1, 1,img_count);
end

fprintf('\tStarting preprocessing: %d images in this stack...\n',img_count)

for cur_img_ind = 1:img_count
    %Find the current input image file
    cur_in_file = find(cur_img_ind <= cumsum(img_count), 1, 'first');
    if cur_in_file > 1
        in_tiff(cur_in_file).setDirectory(cur_img_ind - sum(img_count(1:(cur_in_file-1))));
    else
        in_tiff(cur_in_file).setDirectory(cur_img_ind);
    end
    %Read the image
    cur_img = in_tiff(cur_in_file).read;
    %crop and allign cur img and ref image
    if opts.crop
        %Crop Image
        cur_img = imcrop(cur_img, opts.crop_position);
        if opts.angle <= 90 %if angled to right of yaxis
            cur_img = imrotate(cur_img,-opts.angle);
            ref_img = imrotate(ref_img,-opts.angle);
        else %if angled left of yaxis
            cur_img = imrotate(cur_img,180-opts.angle);
            ref_img = imrotate(ref_img,180-opts.angle);
        end
        cur_img = padarray(cur_img,[60,60]); 
        ref_img = padarray(cur_img,[60 60]);
        ref_img = imcrop(cur_img,[opts.XCropCord,opts.YCropCord,opts.CropHeight,opts.CropWidth]);
        cur_img = imcrop(cur_img,[opts.XCropCord,opts.YCropCord,opts.CropHeight,opts.CropWidth]);
        %regrab the Rfixed image with new size
        Rfixed = imref2d(size(opts.ref_img));
    end
           
    if opts.dsSpatial
        %Spatially downsample 
        cur_img = double(cur_img([1:opts.CropWidth],[1:opts.CropHeight]));
        if opts.ApplyMask
            cur_img = SpatialBin(cur_img,opts.dsFactor,opts.Mask,opts.dsIgnoreNan);
        else
            cur_img = SpatialBin(cur_img,opts.dsFactor);
%             opts.Mask = imresize(opts.Mask,[size(cur_img,1),size(cur_img,2)]);
        end
    else  %Option to mask but not dsSpatial
       if opts.ApplyMask
            cur_img = double(cur_img([1:opts.CropWidth],[1:opts.CropHeight]));
            cur_img(opts.Mask==0) = nan;
       end
    end
    
    stack(:,:,cur_img_ind) = double(cur_img(:,:,1));
    %Write Preprocessed Image to new matstack
    if opts.Verbose
        fprintf('\tFinished preprocessing image %d of %d...\n',cur_img_ind, img_count);
    else
        if mod(cur_img_ind,round(0.02*img_count)) ==0
            fprintf('\t%g%% Complete\n', round(cur_img_ind./img_count*100,2));
        end
    end
end %image loop

%Make DFF;
if opts.makeDFF
[dff,avgproj] = makeDFF(stack, opts);
%save dff
save(opts.DFF_fn,'dff','-V7.3');
end

%Save preprocessed data
save(opts.Stack_fn,'stack','-V7.3');

if opts.motion
    opts.DFF_fn = erase(opts.DFF_fn,'.mat');
    motion_fn = sprintf('%s_registrationinfo.mat',opts.DFF_fn);
    save(motion_fn,'regInfo','-V7.3');
end

fprintf('\tDone with Preprocessing %s...\n',opts.SaveFileBase);
end

