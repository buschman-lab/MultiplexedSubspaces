function [cort_ss,idx] = GetCorticalSubspaces(paired_areas, area_label)
%camden - timeless
%just gets the subspaces between cortical areas
%current only works with 'parent' or 'general' anatomical labels

%get cort areas
cort_label = {'RSPagl','RSPd','RSP','VISp','VISpm','VISa','VISam','VIS','SS','MOs','SSp-bfd','SSp','SSp-un','SSp-n','SSp-m','SSp-tr','SSp-ul','SSp-ll'};

cort_idx = find(ismember(area_label,cort_label)==1);

idx = find(sum(ismember(paired_areas,cort_idx),2)==2);

cort_ss = paired_areas(idx,:);

end %function 

