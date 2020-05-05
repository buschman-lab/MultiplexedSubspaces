
temp = NaN(numel(nanpxs)+size(w,1),size(w,2),size(w,3)); 
goodpxl = 1:size(temp,1);
temp(~ismember(goodpxl,nanpxs),:,:)=w;


figure; 
for i = 1:size(temp,2)
   motif = squeeze(temp(:,i,:));
   motif = reshape(motif,[68,38,size(temp,3)]);
   for j = 1:size(temp,3)
      imagesc(motif(:,:,j),[0 1]);      
      axis equal; axis off
      title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
      pause(0.02);
   end      
end


%%
figure;
% temp = conditionDffMat(data_train',nanpxs,[],[68,38,size(data_train,2)]);
temp = conditionDffMat(temp',nanpxs,[],[68,38,size(temp,2)]);
% temp = conditionDffMat(data_norm',nanpxs,[],[68,38,size(data_norm,2)]);
% temp = conditionDffMat(data_filt,nanpxs,[],[68,38,size(data_filt,1)]);
for j = 1:size(temp,3)
  imagesc(temp(:,:,j),[0 5]);      
  axis equal; axis off
  colormap magma
  title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
  pause();
end      

%%
figure;
% temp = data_norm;
% temp = conditionDffMat(data_norm,nanpxs,[],[68,38,z]);
% temp = dff_frac;
% vals = [prctile(temp(:),2.5),prctile(temp(:),97.5)];
% temp = dff;
for j = 1:size(temp,3)
  imagesc(temp(:,:,j),[0,0.75]);      
  axis equal; axis off
  colormap magma
  title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
  pause(0.001);
end      