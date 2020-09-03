function stack = PreProcess(in_fn,opts)
warning ('off','all'); %Suppress missing tiff metadata warning
%get stack info (faster than tiff counter)
info = imfinfo(in_fn); 
img_count = numel(info); 

%preallocate stack
stack = NaN(ceil(opts.crop_h/opts.spatial_bin_factor),...
    ceil(opts.crop_w/opts.spatial_bin_factor),...
    img_count);

%read tiff
in_tiff = Tiff(in_fn, 'r');

%Loop through each image
for cur_img_ind = 1:img_count
    in_tiff.setDirectory(cur_img_ind);   
    cur_img = in_tiff.read;
    
    %Preprocess image
    cur_img = PreProcessImage(cur_img,opts);

    %Add to stack
    stack(:,:,cur_img_ind) = cur_img; 
    
    %Chatty 
    if opts.verbose
        if mod(cur_img_ind,round(0.01*img_count)) ==0
            fprintf('\t%g%% Complete\n', round(cur_img_ind./img_count*100,2));
        end
    end
    
    
end %image loop
warning ('on','all');


end %function

