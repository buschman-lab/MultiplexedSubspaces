function opts = MaskVasculature(ref_img,opts)

if opts.mask_vasc %mask vasc
    %Indentify vasculature by subtracting a large median filter  (e.g. smoothed image
    %form the base image. This excentuates dark areas that are surrounded by
    %light areas (e.g. vasculature). 
    count=1;
    while 1 %user defined threshold
        if count==1 %use the starting value
            vasc_std = opts.vasc_std;
            neigh = 40; 
            count=count+1;
        else            
            in_vals = inputdlg({'Number of Standard Deviations','Neighboorhood size'},'Manual Threshold',[1,35],{num2str(vasc_std),num2str(neigh)}); %User defined cortical surface
            vasc_std = str2num(in_vals{1});
            neigh = str2num(in_vals{2});
        end
        smoothed = medfilt2(ref_img,[neigh neigh]); %Use large neighboorhood e.g. 1/4 of size
        mask = ref_img-smoothed; 

        %threshold to remove pxs lower that 'x' standard deviations below mean
        %of Vascmask (e.g. the vasculature and other dark blemishes). 
        mask = mask>=(nanmean(nanmean(mask))-(vasc_std*nanstd(nanstd(mask))));

        %close the image to clean up noise
        se = strel('disk',opts.close_disk_size);
        mask = imclose(mask,se); 

        %Add a premade mask to outline the brain.
        if opts.mask_brain_outline 
            temp = load(opts.mask_brain_outline_dir); 
            %Crop the mask to the size of the image. Assumes top left alligned
            temp.brainoutline = temp.brainoutline([1:opts.crop_h],[1:opts.crop_w]);
            mask = mask+temp.brainoutline;
            mask = mask==2;               
        end

        if opts.verbose || opts.manual_mask
            imshowpair(ref_img,ref_img.*mask);
            set(gcf,'units','normalized','position',[0.01 0.04 0.8 0.8]);
            hold on; 
            scatter(opts.bregma(1),opts.bregma(2),'*');
            title(sprintf('Vascmask: nstd = %.2f, neigh = %0.2f',vasc_std,neigh));
        end
        
        temp = questdlg('Threshold look okay?','threshold confirmation');
        if strcmp(temp,'Yes')
            break
        end    
    end %vascular mask threshold
    

    %Manually Nan regions and unmask others 
    if opts.manual_mask
        dlg = questdlg('Would you like to manually mask regions?', ...
        'mask ROIS','Yes','No','No');
        while 1
            switch dlg
                case 'Yes' %mask additional region
                    close all                
                    figure('name','Draw ROI around additional regions to NaN');
                    imshowpair(ref_img,ref_img.*mask);
                    manual_mask = impoly;
                    wait(manual_mask);               
                    manual_mask = manual_mask.createMask();
                    mask(manual_mask)=0; close; 
                    figure('name','Current masked Image')
                    imshowpair(ref_img,ref_img.*mask); 
               case 'No'; break %end dlg
            end
            choice = questdlg(sprintf('Would you like to mask additional regions'),...
                'mask ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
        dlg = questdlg('Would you like to manually UNmask regions?', ...
        'mask ROIS','Yes','No','No');
        while 1 %Manuall UNmask regions 
            switch dlg
                case 'Yes' %mask additional region
                    close all
                    figure('name','Draw ROI around additional regions to NaN');
                    imshowpair(ref_img,ref_img.*mask);
                    manual_mask = impoly;
                    wait(manual_mask);               
                    manual_mask = manual_mask.createMask(); %his is 1 in the desired ROI          
                    mask(manual_mask)=1; close
                    figure('name','Current masked Image')
                    imshowpair(ref_img,ref_img.*mask);
               case 'No'; break
            end
            choice = questdlg(sprintf('Would you like to unmask additional regions'),...
                'mask ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
    end

    if opts.verbose || opts.manual_mask
        close all
        imshowpair(ref_img,ref_img.*mask);
        hold on; 
        scatter(opts.bregma(1),opts.bregma(2),'*');
        title('FINAL cropped Img with Vascmask overlay');
            choice = questdlg(sprintf('Would you like to flag this brain for imperfections?'),...
            'Flag brain',...
            'Yes','No','No');
       switch choice
           case 'Yes'
               opts.flagged_brain = 1;
           case 'No' 
               opts.flagged_brain = 0; 
       end
    elseif opts.maskVasculature==0 && opts.maskBrainOutline == 1 %only apply exterior mask
        mask = brainoutline;
        if opts.verbose
            close all
            imshowpair(ref_img,ref_img.*mask);
            hold on; 
            scatter(opts.bregma(1),opts.bregma(2),'*');
            title('FINAL cropped Img with masked Brain Overlay');
        end        
    else %apply no mask
        mask = ones(size(ref_img));
    end
    opts.mask = mask; 
    opts.masked_ref_img = ref_img.*mask; 

end %function end
    