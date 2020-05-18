function [tform,output_size]= RegisterImages(fixed,moving,type)
%Camden MacDowell 2020
%toolbox of different registration methods

%first preprocess image
fixed = (imadjust(uint16(fixed)));
moving = (imadjust(uint16(moving)));

switch type
    case 'auto'
        %First try auto registration
        threshold = 1000;
        while 1   %iteratively decrease the threshold, until enough points for matching (4)
            %detect features
            pts_original  = detectSURFFeatures(fixed,'MetricThreshold',threshold,'NumOctaves',4);
            pts_new = detectSURFFeatures(moving,'MetricThreshold',threshold);
            [features_original,  valid_pts_original]  = extractFeatures(fixed,  pts_original);
            [features_new, valid_pts_new ] = extractFeatures(moving, pts_new);

            %match features
            indexPairs = matchFeatures(features_original, features_new);
            matched_original = valid_pts_original(indexPairs(:,1));
            matched_new = valid_pts_new (indexPairs(:,2));

            if size(matched_new,1) < 2
                threshold = threshold - 50;
                if threshold <= 0
                    error('Image Registration Failed: Threshold set to zero');
                end
            else       
                %get transformation
                [tform, ~, ~] = estimateGeometricTransform(matched_new , matched_original, 'similarity');
                output_size = imref2d(size(fixed));
                % registered = imwarp(distorted,tform,'OutputView',output_size);
                break
            end
        end        
    case 'manual' %use manually selected fudicial points to register images
        [cp_moving,cp_fixed] = cpselect(moving,fixed,'Wait',true);
        %fine tune with cross correlation
        cp_moving_adjusted = cpcorr(cp_moving,cp_fixed,moving,fixed);        
        tform = fitgeotrans(cp_moving_adjusted,cp_fixed,'similarity');
        output_size = imref2d(size(fixed)); %relate intrinsic and world coordinates
%         registered = imwarp(moving,tform,'OutputView',output_size); 
%         figure
%         imshowpair(registered,fixed,'blend')        
end
 
%     %show user, validation images
%     figure;
%     showMatchedFeatures(original,new,matched_original,matched_new);
%     title('Putatively matched points (including outliers)');
%     %Visualize outcome of registration
%     registered_img = PreProcessImage(new,opts,{'spatial_bin_factor',1});
%     figure; hold on; imagesc(registered_img); colormap('gray'); axis off
%     set(gca,'ydir','reverse','xlim',[1,size(registered_img,1)],'ylim',[1,size(registered_img,2)]);
%     plot(prepro_log.bregma(1),opts.bregma(2),'x','markersize',10,'color','r','LineWidth',3);
%     
%     %ask if okay, otherwise do manual registration
%     while 1
%         dlg = questdlg('Valid Regiration?', ...
%         'Registration','Yes','No','No');
%         switch dlg
%             case 'No' %try manual registration
%                 close all                
% 
%            case 'Yes'; break %end dlg
%         end
%     end

% 
% %select up to 10 input points
%         figure('units','normalized','position',[0 0 1 1]); hold on;          
%         ax_fixed = subplot(2,1,1); imshow(fixed); axis off; hold on; 
%         ax_moving = subplot(2,1,2); imshow(moving); axis off; hold on; 
%         title('select points on fixed image right click when done','fontsize',10);
%         
%         %work on the fixed image
%         axes(ax_fixed);
%         COUNT = 1; cp_fixed = ones(10,2); 
%         while 1   % read ginputs until a mouse right-button occurs
%             [cp_fixed(COUNT,1),cp_fixed(COUNT,2),button] = ginput(1);
%             if button <=1
%                 plot(cp_fixed(COUNT,1),cp_fixed(COUNT,2),'x','color','r','linewidth',2,'MarkerSize',10);
%                 text(cp_fixed(COUNT,1)+2,cp_fixed(COUNT,2)+2,num2str(COUNT),'color','r','FontSize',12,'Fontweight','normal');
%                 COUNT = COUNT+1;
%             else %remove the last point (due to clicking off the screen
%                 cp_fixed(COUNT:end,:) = [];
%                 break
%             end                
%         end
%         
%         %moving image
%         axes(ax_moving);
%         %select points on target image  
%         title('Now select same points on moving image','fontsize',10);
%         cp_moving = ones(size(cp_fixed,1),2);
%         for i = 1:size(cp_fixed,1)
%             [cp_moving(i,1),cp_moving(i,2)] = ginput(1);
%             plot(cp_moving(i,1),cp_moving(i,2),'x','color','g','linewidth',2,'MarkerSize',10);
%             text(cp_moving(i,1)+2,cp_moving(i,2)+2,num2str(i),'color','g','FontSize',12,'Fontweight','normal');
%         end     
%         close;
%         
%         %fine tune with cross correlation
%         cp_moving_adjusted = cpcorr(cp_moving,cp_fixed,moving,fixed);
%         
%         %confirm
%         figure(); hold on;          
%         imshow(moving); axis off; hold on; 
%         title('Auto-adjusted points using cross correlation','fontsize',10);
%         plot(cp_moving_adjusted(:,1),cp_moving_adjusted(:,2),'x','color','c','MarkerSize',10,'linewidth',2)      
%         plot(cp_moving(:,1),cp_moving(:,2),'x','color','g','MarkerSize',10,'linewidth',2)     


