function [bregma] = setBregma(ref_img)

% user select bregma
figure('name','select bregma','units','normalized','outerposition',[0 0 1 1]);
imagesc(ref_img); colormap gray; axis image
bregma = round(ginput(1));
close

end