function PlotProbeAnatomy(axis_handle, neu_area, xpos, type,direction,offset)
%Camden MacDowell - timeless
%plot a multicolor line at a given x position with colored text labels next
%to it (clipping off)

%EDIT: needs to be ordered in the desired order prior to feeding in to function
% neu_depth = clust_info.depth;
% %ordered index by surface to tip (increasign)
% [~,idx] = sort(neu_depth,'ascend');
% neu_area = neu_area(idx);

%axes(axis_handle)
switch type
    case 'parent'
        label = [neu_area.parent_label];
    case 'detail'
        label = [neu_area.detailed_label];
    otherwise
        error('unknown type');
end

axes(axis_handle); hold on; 
unique_label = unique(label); 
c = distinguishable_colors(numel(unique_label));
for i = 1:numel(unique_label)
    yvals = find(strcmp(label,unique_label{i})==1);    
    %if multiple words, split to new lines
    temp = regexp(unique_label{i},' ','split');      
    if direction==1
        plot([xpos,xpos],[yvals(1),yvals(end)],'color',c(i,:),'linewidth',2)   
        text(xpos-offset,[yvals(1)+(yvals(end)-yvals(1))/2],sprintf('%s\n',temp{:}),'color',c(i,:),'Rotation',0,'HorizontalAlignment','right')
    else
        plot([yvals(1),yvals(end)],[xpos,xpos],'color',c(i,:),'linewidth',2)   
        text([yvals(1)+(yvals(end)-yvals(1))/2],xpos-offset,sprintf('%s\n',temp{:}),'color',c(i,:),'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','top')
    end
end

set(gca,'Clipping','off');

end



