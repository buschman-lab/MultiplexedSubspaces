function [angle] = setMidlineAngle(meanproj)
% user-set bregma coordinates + midline for rotation 
% Modified from Lucas Pinto 2016 by Camden MacDowell 2018

% user select two points to fit midline
figure('name','select two points for midline','units','normalized','outerposition',[0 0 1 1]);
imagesc(meanproj); colormap gray; axis image
midline = imline; %drawline functionality in 2018b isn't quite there yet so don't use
position = wait(midline);
close
% as in y = ax + b
a = (position(1,2)-position(2,2))/(position(1,1)-position(2,1));
b = position(1,2) - a*position(1,1);

%Rotation of the brain with respect to y axis
angle = 90-(atan(a) * (180 / pi));

%save bregma
end