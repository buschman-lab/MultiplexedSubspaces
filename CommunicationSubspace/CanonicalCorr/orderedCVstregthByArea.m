function orderedCVstregthByArea(data)
%Camden - timeless
%Creates a plot per area showing the strength of subspaces per motif


%Area list
area_label = data(1).area_label;
paired_areas = data(1).paired_areas;
[paired_areas,cort_idx] = GetCorticalSubspaces(paired_areas, area_label);


%get the strength for each motif
[~,x] = cvStrengthMap(data,'r_all_sum');
x(:,2)=[];
%loop through each subspace
for i = 1:size(paired_areas,1)   
   [temp,idx] = sort(x(i,:),'ascend');
   temp(isnan(temp))=0;
   figure; hold on; 
   title(sprintf('%s and %s',area_label{paired_areas(i,1)},area_label{paired_areas(i,2)}),'fontweight','normal')
   plot(1:numel(temp),temp,'marker','x','linestyle','none')
   idx = arrayfun(@(x) num2str(x), idx,'UniformOutput',0);
   text(1:numel(temp),temp,idx)
end

end



