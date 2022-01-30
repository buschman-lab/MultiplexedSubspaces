function [area_val, area_label] = CleanUpAreas(area_val, area_label, min_num)
%Camden - timeless
%Cleans up any erraneous spikes in fiber tracks or unlabeled regions
%pass in min_num to remove areas with less than min_num neurons

if nargin <3; min_num = 0; end

%remove fiber tracts and unlabeled areas
bad_label = {'cc','na','fiber tracts'}; 
bad_idx = cellfun(@(x) find(strcmp(area_label,x)==1),bad_label,'UniformOutput',0);
bad_idx(cellfun(@ isempty, bad_idx))=[];
bad_idx = [bad_idx{:}]; 

area_val(bad_idx) = [];
area_label(bad_idx) = [];

%remove with too few
n = cellfun(@(x) size(x,1), area_val);
area_val(n<min_num)=[];
area_label(n<min_num)=[];

end