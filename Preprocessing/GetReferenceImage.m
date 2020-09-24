function ref_img = GetReferenceImage(in_fn,fixed_image)
%FixedImage: 0 if first; 1 if middle, 2 if mean image (slow)

%if available use the TIFFStack Module (way faster)
%https://github.com/DylanMuir/TIFFStack
if exist('TIFFStack','class')==8    
    warning ('off','all'); %turn off because of header warnings
    tsStack = TIFFStack(in_fn);    
    
    [~,~,nZ] = size(tsStack);
    %Grab the first image
    switch fixed_image        
        case 'middle'
            ref_img = double(tsStack(:,:,floor(nZ/2)));
        case 'mean'
            %decimate because no need to take whole average
            if nZ>100
                ref_img = nanmean(tsStack(:,:,1:20:end),3);
            else
                ref_img = nanmean(tsStack,3);
            end
        case 'max'
            %decimate because no need to take whole average
            if nZ>100
                ref_img = nanmax(tsStack(:,:,1:20:end),[],3);
            else
                ref_img = nanmax(tsStack,[],3);
            end       
        case 'variance'
            %decimate because no need to take whole average
            if nZ>100
                ref_img = nanvar(single(tsStack(:,:,1:20:end)),[],3);
            else
                ref_img = nanvar(single(tsStack),[],3);
            end       
        otherwise  %use first image
            warning ('on','all');
            warning('unknown image selection, using first image')
    end
            
    warning ('on','all');
    
 


%%%%LEGECY%%%%%

else %LEGaCY CODE: If not @TIFFStack Object oriented module
    %Surpress warning regarding tag structure of tiff. 
    warning ('on','all');
    warning('Using Legacy Code To Load Data (omg... so sloooow)')
    warning ('off','all');
    in_tiff = Tiff(in_fn, 'r');

    %Count how many images in each file
    in_tiff.setDirectory(1);
    while ~in_tiff.lastDirectory
        in_tiff.nextDirectory;
    end
    img_count = in_tiff.currentDirectory;

    %stop from going on for a long time
    if img_count >500
       img_count = 500;
    end

    %Grab the middle image
     switch fixed_image 
         case 'middle'
            img_idx = floor(sum(img_count)/2);

            fprintf('\nUsing middle image as ref image (%d of %d)\n', img_idx, img_count);
            cur_in_file = find(img_idx <= cumsum(img_count), 1, 'first');
            if cur_in_file > 1
                in_tiff(cur_in_file).setDirectory(img_idx - sum(img_count(1:(cur_in_file-1))));
            else
                in_tiff(cur_in_file).setDirectory(img_idx);
            end
            ref_img = in_tiff(cur_in_file).read;
        
         case 'mean'
            fprintf('\nCalculating average image of first % to use as fixed image...\n');
            fixed_img = zeros(in_tiff(1).getTag('ImageLength'),in_tiff(1).getTag('ImageWidth'));
            stack = uint16(zeros(size(fixed_img,1),size(fixed_img,2),img_count));
            for cur_img_ind = 1:2:img_count
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
             
        otherwise
            img_idx = 1;
            fprintf('\nUsing first image as ref image');
            cur_in_file = find(img_idx <= cumsum(img_count), 1, 'first');
            if cur_in_file > 1
                in_tiff(cur_in_file).setDirectory(img_idx - sum(img_count(1:(cur_in_file-1))));
            else
                in_tiff(cur_in_file).setDirectory(img_idx);
            end
            ref_img = in_tiff(cur_in_file).read;
    end
    %turn on warnings again
    warning ('on','all');

end