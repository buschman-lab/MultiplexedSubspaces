function [data, labels, Perc_Filtered] = parse_dlc(raw_data,parts_list,reference_part,epsilon)
%dlc stores data as parts and then x,y,likihood. 
%and trims at row 4 

%remove all columns with likelihood values; 
idx = regexp(raw_data(3,:),'likelihood','match');
idx = cellfun(@(x) size(x,1),idx,'UniformOutput',0);
raw_data(:,[idx{:}]==1)=[];

if ~isempty(reference_part)
    idx = regexp(raw_data(2,:),reference_part,'match');
    idx = cellfun(@(x) size(x,1),idx,'UniformOutput',0);
    reference = raw_data(4:end,[idx{:}]==1);
    reference = cell2mat(reference);
end

%loop through parts list and get the columns for each 
idx = cellfun(@(x) regexp(raw_data(2,:),x),parts_list,'UniformOutput',0);
for i = 1:numel(idx)
    temp = cellfun(@(x) ~isempty(x), idx{i}, 'UniformOutput',0);
    idx{i} = [temp{:}];    
end
idx = sum(cat(1,idx{:}),1);

%gut checks
if any(idx>1)
    error('you have repeats in the parts list');
end

if sum(idx)<2*numel(parts_list)
    warning('not all parts found');
end

data = raw_data(4:end,idx==1);    
data = cell2mat(data);

%optional filter 
if epsilon
    %filter the reference
    if ~isempty(reference_part)
        [reference, ~] = filter_dlc(reference,epsilon);
    end
    %filter the data
    [data, change_idx] = filter_dlc(data,epsilon);
    Perc_Filtered = nansum(change_idx)/numel(change_idx);
else
    Perc_Filtered = NaN;
end

%perform function using reference channel
if ~isempty(reference_part)
    data(:,1:2:end) = data(:,1:2:end)-reference(:,1);
    data(:,2:2:end) = data(:,2:2:end)-reference(:,2);
end

labels = repmat(parts_list,2,1);
labels = labels(:);

end




function [data, change_idx] = filter_dlc(data,epsilon)
    change_idx = NaN(size(data));
    for col = 1:size(data,2)
        for row = 2:size(data,1)
           temp = abs(diff([data(row,col),data(row-1,col)]));
           if temp>=epsilon
               data(row,col)=data(row-1,col); %replace pixels that move too much with the previous location
               change_idx(row,col) = 1;
           end
        end
    end    
end















