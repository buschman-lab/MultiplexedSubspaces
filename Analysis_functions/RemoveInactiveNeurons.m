function [area_val, inactive_idx] = RemoveInactiveNeurons(area_val, min_fr)
%Camden - timeless
%Cleans up any erraneous spikes in fiber tracks or unlabeled regions
%pass in min_num to remove areas with less than min_num neurons

if nargin <2; min_fr = 0.5/15; end

inactive_idx = cellfun(@(x) nanmean(x,[2,3])<=min_fr, area_val,'UniformOutput',0);

area_val = cellfun(@(x,y) x(y==0,:,:),area_val,inactive_idx,'UniformOutput',0);


end