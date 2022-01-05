function PlotProbeAnatomy(axis_handle, neu_area, xpos, type,direction, offset, c)
%Camden MacDowell - timeless
%plot a multicolor line at a given x position with colored text labels next
%to it (clipping off)

%EDIT: needs to be ordered in the desired order prior to feeding in to function
% neu_depth = clust_info.depth;
% %ordered index by surface to tip (increasign)
% [~,idx] = sort(neu_depth,'ascend');
% neu_area = neu_area(idx);

if nargin <7 %legacy contingency
    c = [];
end

%axes(axis_handle)
idx = find(arrayfun(@(n) isempty(neu_area(n).parent_label), 1:size(neu_area,2),'UniformOutput',1)==1); %fill in unlabeled regions
for j = idx        
    neu_area(j).detailed_label = {''};
    neu_area(j).parent_label = {''};
end
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
if isempty(c)
   c = distinguishable_colors(numel(unique_label));
end

%check for and split, non contiguous trajectories into the same 
%region as different entries
yvals_all = {};
rep_label = {};
c_all = {};
COUNT = 1;
for i = 1:numel(unique_label)
   temp_yvals = find(strcmp(label,unique_label{i})==1);
   split_idx = find([0 diff(temp_yvals)>1]==1);
   if ~isempty(split_idx)
       for j = 1:numel(split_idx)
           if j == 1
               yvals_all{COUNT} = temp_yvals(1:split_idx(j)-1);
           else
               yvals_all{COUNT} = temp_yvals(split_idx(j-1):split_idx(j)-1);
           end       
           rep_label{COUNT} = unique_label{i};
           c_all{COUNT} = c(i,:);
           COUNT = COUNT+1;
       end
       yvals_all{COUNT} = temp_yvals(split_idx(j):end); %add the final section
       rep_label{COUNT} = unique_label{i};
       c_all{COUNT} = c(i,:);
       COUNT = COUNT+1;       
   else
       yvals_all{COUNT} = temp_yvals;
       rep_label{COUNT} = unique_label{i};
       c_all{COUNT} = c(i,:);
       COUNT = COUNT+1;
   end    
end
c_all = cat(1,c_all{:});

for i = 1:numel(rep_label)   
    yvals = yvals_all{i};
    %if multiple words, split to new lines
%     if ~isempty(regexp(rep_label{i},' ','match'))
    temp = regexp(rep_label{i},' ','split');  
%     elseif ~isempty(regexp(rep_label{i},'-','match'))
%         temp = regexp(rep_label{i},'-','split');          
%     else
%         temp = rep_label(i);  
%     end
    if direction==1
        plot([xpos,xpos],[yvals(1)-0.5,yvals(end)+0.5],'color',c_all(i,:),'linewidth',2)   
        text(xpos-offset,[yvals(1)+(yvals(end)-yvals(1))/2],sprintf('%s\n',temp{:}),'color',c_all(i,:),'Rotation',0,'HorizontalAlignment','right')
    else
        plot([yvals(1)-0.5,yvals(end)+0.5],[xpos,xpos],'color',c_all(i,:),'linewidth',2)   
        text([yvals(1)+(yvals(end)-yvals(1))/2],xpos-offset,sprintf('%s\n',temp{:}),'color',c_all(i,:),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle')
    end
end

set(gca,'Clipping','off');

end



