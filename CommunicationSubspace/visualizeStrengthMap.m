function fh = visualizeStrengthMap(cv_map,data,c_val)
%Camden - timeless

if nargin <3; c_val = []; end
area_label = data(1).area_label;
paired_areas = data(1).paired_areas;

%if not in matrix form
if size(cv_map,1)~=size(cv_map,2)
    cv_map = SubspaceStatMap(cv_map,paired_areas);
end

figure('units','normalized','position',[0 0 1 1]); hold on; 
[num_rows, num_col]=numSubplot(size(cv_map,3),2);
for i = 1:size(cv_map,3)
   subplot(num_rows,num_col,i); hold on; 
   if ~isempty(c_val)
       imagesc(cv_map(:,:,i),c_val); 
   else
       imagesc(cv_map(:,:,i)); 
   end
   colormap(flipud(gray))
   set(gca,'xtick',1:numel(area_label),'XTickLabel',area_label,'XTickLabelRotation',45);
   set(gca,'ytick',1:numel(area_label),'YTickLabel',area_label);
   colorbar
   axis square
   box off
   set(gca,'xlim',[0.5 size(cv_map,2)+0.5]);
   set(gca,'ylim',[0.5 size(cv_map,2)+0.5])
   title(sprintf('motif %d',i)); 
end

fh = gcf; 

end
