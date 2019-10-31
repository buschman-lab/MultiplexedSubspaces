function roi = SelectROI(img,label,other_roi)

if nargin<3; other_roi = []; end %other_roi is a cell array of rois to also show on the plot
%show image
figure('name',sprintf('Select %s', label))
imagesc(img); colormap gray

if ~isempty(other_roi)
    cellfun(@(x) rectangle('Position',x.position,'EdgeColor',x.color,'FaceColor',x.color),other_roi,'UniformOutput',0);
end

%Set to correct aspect ratio and shift to ~center of 1080p screen
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1)-200, pos(2)-200, size(img,2), size(img,1)]);
hold on

%cropping rectangle
hROI = drawrectangle(gca);
roi = CustomWait(hROI);

%add label to structure
roi.label = label; 

close

end %function end



function roi = CustomWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
roi.position = round(hROI.Position,0);
roi.color = hROI.Color;

end



function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end