function MakeGif(data,filename,delay_time,marker_index)
if nargin <4; marker_index = 0; end
    

v = VideoWriter([filename '.avi'],'Motion JPEG AVI');
v.FrameRate = 60;
open(v)
for i = 1:size(data,3)
    imagesc(data(:,:,i));
    axis off
    colormap gray    
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame)
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File 
%     if i == 1 
%         imwrite(frame,cm,filename,'gif', 'Loopcount',1,'DelayTime',delay_time); 
%     else 
%         imwrite(frame,cm,filename,'gif','WriteMode','append','DelayTime',delay_time); 
%     end 
end
close(v)

end
