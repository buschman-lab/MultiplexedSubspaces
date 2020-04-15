function cur_img = PreProcessImage(cur_img,opts,myvarargin)
if nargin <3; myvarargin = []; end
opts = ParseOptionalInputs(opts,myvarargin);

%Register image (optional)
if ~isempty(opts.tform)
   cur_img = imwarp(cur_img,opts.tform,'OutputView',opts.output_size);
end

%Crop and Allign Image
cur_img = imcrop(cur_img, opts.crop_position);
if opts.angle <= 90 %if angled to right of yaxis
    cur_img = imrotate(cur_img,-opts.angle);
else %if angled left of yaxis
    cur_img = imrotate(cur_img,180-opts.angle);
end
cur_img = padarray(cur_img,[60,60]); 
cur_img = imcrop(cur_img,[opts.crop_cord(1),opts.crop_cord(2),...
    opts.crop_w,opts.crop_h]);

%confirm correct dimensions since imcrop is inconsistent (see imcrop)
cur_img = double(cur_img([1:opts.crop_h],[1:opts.crop_w])); 

%Mask and spatial bin
cur_img = SpatialBin(cur_img,opts.spatial_bin_factor,opts.mask,1);

end