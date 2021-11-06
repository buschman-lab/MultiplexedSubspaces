function [activity_in,activity_contra]=CompareROIVariance(dff,probe_coords,type)
%Camden MacDowell - timeless
%plots the histograms of pixels values in each craniotomy
%type can be pixelwise variance ('pxlvar') or all values over time ('all')

R_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
R_hemi(:,37:end,:) = fliplr(dff(:,3:34,:));
R_hemi = reshape(R_hemi,[size(R_hemi,1)^2,size(R_hemi,3)]);
L_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
L_hemi(:,37:end,:) = dff(:,37:end,:);
L_hemi = reshape(L_hemi,[size(L_hemi,1)^2,size(L_hemi,3)]);

mask = zeros(size(dff,1),size(dff,2));

within_mask = cell(1,4);
for i = 1:numel(probe_coords)
    temp = probe_coords{i};
    within_mask{i} = mask;
    within_mask{i}(temp(1,2)+[-2:2],temp(1,1)+[-2:2])=1; %remember, x and y need to be flipped from what you would expect
    within_mask{i} = reshape(within_mask{i},[size(mask,1)^2,size(mask,3)]); %flatten
end

switch type 
    case 'all'
        activity_in = cellfun(@(x) nanmean(L_hemi(x==1,:))', within_mask,'UniformOutput',0);
        activity_in = [activity_in{:}];
        activity_contra = cellfun(@(x) nanmean(R_hemi(x==1,:))', within_mask,'UniformOutput',0);
        activity_contra = [activity_contra{:}];
    case 'pxlstd' %START HERE: 
        activity_in = cellfun(@(x) nanstd(L_hemi(x==1,:),[],2), within_mask,'UniformOutput',0);
        activity_in = [activity_in{:}];
        activity_contra = cellfun(@(x) nanstd(R_hemi(x==1,:),[],2), within_mask,'UniformOutput',0);
        activity_contra = [activity_contra{:}];
    case 'pxlvar'
        activity_in = cellfun(@(x) nanvar(L_hemi(x==1,:),[],2), within_mask,'UniformOutput',0);
        activity_in = [activity_in{:}];
        activity_contra = cellfun(@(x) nanvar(R_hemi(x==1,:),[],2), within_mask,'UniformOutput',0);
        activity_contra = [activity_contra{:}];     
    otherwise
        error('unknown type');
end


end

