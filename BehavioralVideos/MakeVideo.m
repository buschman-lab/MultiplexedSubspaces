function MakeVideo(data,filename,marker_index)
if nargin <3; marker_index = 0; end
    

v = VideoWriter([filename '.avi'],'Motion JPEG AVI');
v.FrameRate = 60;
open(v)
for i = 1:size(data,3)
    cla
    imagesc(data(:,:,i));
    axis off
    colormap gray    
    if ismember(marker_index,i)
       hold on; scatter(20,20,50,'r');
    end
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v)

end
