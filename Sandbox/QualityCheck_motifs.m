
temp = NaN(numel(nanpxs)+size(w,1),size(w,2),size(w,3)); 
goodpxl = 1:size(temp,1);
temp(~ismember(goodpxl,nanpxs),:,:)=w;


figure; 
active = find(nanmean(squeeze(nanmean(w,3)),1)>0);
for i = active
   motif = squeeze(temp(:,i,:));
   motif = reshape(motif,[68,38,size(temp,3)]);
   for j = 1:size(temp,3)
      imagesc(motif(:,:,j),[0 1.5]); colorbar     
      axis equal; axis off
      title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
      pause();
   end      
end

%%
temp = NaN(numel(nanpxs)+size(basis_motifs,1),size(basis_motifs,2),size(basis_motifs,3)); 
goodpxl = 1:size(temp,1);
temp(~ismember(goodpxl,nanpxs),:,:)=basis_motifs;

figure; 
active = find(nanmean(squeeze(nanmean(basis_motifs,3)),1)>0);
for i = active
   motif = squeeze(temp(:,i,:));
   motif = reshape(motif,[68,38,size(temp,3)]);
   for j = 1:size(temp,3)
      imagesc(motif(:,:,j),[0 1]); colorbar     
      axis equal; axis off
      title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
      pause();
   end      
end


%%
figure;
% temp = conditionDffMat(data_train',nanpxs,[],[68,38,size(data_train,2)]);
% temp = conditionDffMat(temp',nanpxs,[],[68,38,size(temp,2)]);
% temp = conditionDffMat(data_norm',nanpxs,[],[68,38,size(data_norm,2)]);
% temp = conditionDffMat(data_filt,nanpxs,[],[68,38,size(data_filt,1)]);
temp = conditionDffMat(data,nanpxs,[],[68,38,size(data,1)]);
for j = 1:size(temp,3)
  imagesc(temp(:,:,j),[0 0.75]);      
  axis equal; axis off
  colormap magma
  title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
  pause(0.001);
end      

%%
figure;
% temp = data_norm;
% temp = conditionDffMat(data_norm,nanpxs,[],[68,38,z]);
% temp = dff_frac;
% vals = [prctile(temp(:),2.5),prctile(temp(:),97.5)];
% temp = dff;
for j = 7000:size(temp,3)
  imagesc(temp(:,:,j),[0,0.5]);      
  axis equal; axis off
  colormap magma
  title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
  pause(0.001);
end      