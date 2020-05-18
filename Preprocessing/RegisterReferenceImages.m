function opts = RegisterReferenceImages(fixed,moving,opts)
%try different types of image registration to see which works. loop through
%user attempts at manual registration if automated fails

%%
%first try to autoregistration (fast)
try
    [opts.tform,opts.output_size] = RegisterImages(fixed,moving,'auto'); 
catch
    [opts.tform,opts.output_size] = RegisterImages(fixed,moving,'manual'); 
end
ValidationImages(moving,opts);
%%
%ask if okay, otherwise do manual registration
while 1
    dlg = questdlg('Valid Regiration?', ...
    'Registration','Yes','No','No');
    switch dlg
        case 'No' %try manual registration
            close all                
            [opts.tform,opts.output_size] = RegisterImages(fixed,moving,'manual'); 
            ValidationImages(moving,opts);
       case 'Yes'; break %end dlg
    end
end

end


function ValidationImages(moving,opts)
    %Visualize outcome of registration
    registered_img = PreProcessImage(moving,opts,{'spatial_bin_factor',1});    
    figure; hold on; imagesc(registered_img); colormap('gray'); axis off; axis equal
    set(gca,'ydir','reverse','ylim',[1,size(registered_img,1)],'xlim',[1,size(registered_img,2)]);
    plot(opts.bregma(1),opts.bregma(2),'x','markersize',10,'color','r','LineWidth',3);
end   
 