function [area_val, area_label] = ParseByArea(x,neu_area,type)
%Camden MacDowell - timeless
%takes any input vector of equal size to number of units and parses by the
%unique anatomical areas in neu_area

switch type
    case 'detail'
        temp = cat(1,neu_area(:).detailed_label);
    case 'parent'
        temp = cat(1,neu_area(:).parent_label);
    otherwise
        error('unknown level of anatomical detail');
end

%split by unique areas
area_label = unique(temp);
area_val = cellfun(@(y) x(strcmp(temp,y),:,:,:), area_label,'UniformOutput',0); %supports up to 4D


end %function