function [data, Perc_Filtered] = parse_dlc(raw_data,bp)
if nargin <2; bp = behavioral_params; end

%dlc stores data as parts and then x,y,likihood. 
%and trims at row 4 

%remove all columns with likelihood values; 
idx = regexp(raw_data(3,:),'likelihood','match');
idx = cellfun(@(x) size(x,1),idx,'UniformOutput',0);
raw_data(:,[idx{:}]==1)=[];

if  ~isemtpy(bp.dlc_reference_part)
    idx = regexp(raw_data(2,:),reference_part,'match');
    idx = cellfun(@(x) size(x,1),idx,'UniformOutput',0);
    reference = raw_data(4:end,[idx{:}]==1);
    reference = cell2mat(reference);
end

%loop through parts list and get the columns for each 
idx = cellfun(@(x) regexp(raw_data(2,:),x),bp.dlc_parts_list,'UniformOutput',0);
for i = 1:numel(idx)
    temp = cellfun(@(x) ~isempty(x), idx{i}, 'UniformOutput',0);
    idx{i} = [temp{:}];    
end
idx = sum(cat(1,idx{:}),1);

%gut checks
if any(idx>1)
    error('you have repeats in the parts list');
end

if sum(idx)<2*numel(bp.dlc_parts_list)
    warning('not all parts found');
end

data = raw_data(4:end,idx==1);    
data = cell2mat(data);

%optional filter 
if bp.dlc_epsilon ~=0
    %filter the reference
    [reference, ~] = filter_dlc(reference,bp.dlc_epsilon);
    %filter the data
    [data, change_idx] = filter_dlc(data,bp.dlc_epsilon);
    Perc_Filtered = nansum(change_idx)/numel(change_idx);
else
    Perc_Filtered = NaN;
end

%subtract reference channel
data(:,1:2:end) = data(:,1:2:end)-reference(:,1);
data(:,2:2:end) = data(:,2:2:end)-reference(:,2);

end




function [data, change_idx] = filter_dlc(data,epsilon)
    change_idx = NaN(size(data));
    temp = abs(diff(data,1));
    for col = 1:size(temp,2)
        for row = 2:size(temp,1)
           if temp(row,col)>=epsilon
               data(row,col)=data(row-1,col); %replace pixels that move too much with the previous location
               change_idx(row,col) = 1;
           end
        end
    end    
end















