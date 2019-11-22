function stack = PreProcess(in_fn,opts)
warning ('off','all'); %Suppress missing tiff metadata warning
%get stack info (faster than tiff counter)
info = imfinfo(in_fn); 
img_count = numel(info); 

%preallocate stack
stack = NaN(floor(opts.crop_w/opts.spatial_bin_factor),...
    floor(opts.crop_h/opts.spatial_bin_factor),...
    img_count);

%read tiff
in_tiff = Tiff(in_fn, 'r');

%Loop through each image
for cur_img_ind = 1:img_count
    in_tiff.setDirectory(cur_img_ind);   
    cur_img = in_tiff.read;
    
    %Crop and Allign Image
    cur_img = imcrop(cur_img, opts.crop_position);
    if opts.angle <= 90 %if angled to right of yaxis
        cur_img = imrotate(cur_img,-opts.angle);
    else %if angled left of yaxis
        cur_img = imrotate(cur_img,180-opts.angle);
    end
    cur_img = padarray(cur_img,[60,60]); 
    cur_img = imcrop(cur_img,[opts.crop_cord(1),opts.crop_cord(2),...
        opts.crop_h,opts.crop_w]);
    %confirm correct dimensions since imcrop is inconsistent (see imcrop)
    cur_img = double(cur_img([1:opts.crop_w],[1:opts.crop_h])); 
        
    %Mask and spatial bin
    cur_img = SpatialBin(cur_img,opts.spatial_bin_factor,opts.mask,1);

    %Add to stack
    stack(:,:,cur_img_ind) = cur_img; 
    
    %Chatty 
    if opts.verbose
        if mod(cur_img_ind,round(0.02*img_count)) ==0
            fprintf('\t%g%% Complete\n', round(cur_img_ind./img_count*100,2));
        end
    end
    
    
end %image loop
warning ('on','all');


end %function

