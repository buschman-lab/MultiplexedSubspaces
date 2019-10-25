function stack = PreProcess(in_fn,opts)
%get stack info (faster than tiff counter)
info = imfinfo(in_fn); 
img_count = numel(info); 

%preallocate stack
stack = NaN(floor(opts.crop_w/opts.spatial_bin_factor),...
    floor(opts.crop_h/opts.spatial_bin_factor),...
    img_count);

%read tiff
in_tiff = Tiff(in_fn, 'r');

%Suppress missing tiff metadata warning
warning ('off','all');
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


% 
% 
% %Make DFF;
% if opts.makeDFF
% [dff,avgproj] = makeDFF(stack, opts);
% %save dff
% save(opts.DFF_fn,'dff','-V7.3');
% end
% 
% %Save preprocessed data
% save(opts.Stack_fn,'stack','-V7.3');
% 
% if opts.motion
%     opts.DFF_fn = erase(opts.DFF_fn,'.mat');
%     motion_fn = sprintf('%s_registrationinfo.mat',opts.DFF_fn);
%     save(motion_fn,'regInfo','-V7.3');
% end
% 
% fprintf('\tDone with Preprocessing %s...\n',opts.SaveFileBase);
% % if opts.masc_mask == 1
%     stack = zeros(size(opts.Mask,1),size(opts.Mask,2));
% else
%     if opts.crop
%         stack = zeros(floor(opts.crop_position(3)),floor(opts.crop_position(4)));
%     else
%         stack = zeros(floor(size(ref_img,1)),floor(size(ref_img,2)));
%     end
% end
% %Preallocate along the 3rd dimension
% if opts.dsSpatial
%     stack = SpatialBin(stack,opts.dsFactor);
% else
%     stack = repmat(stack, 1, 1,img_count);
% end
