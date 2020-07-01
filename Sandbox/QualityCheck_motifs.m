% w = W;
temp = NaN(numel(nanpxs)+size(w,1),size(w,2),size(w,3)); 
goodpxl = 1:size(temp,1);
temp(~ismember(goodpxl,nanpxs),:,:)=w;


figure; 
active = find(nanmean(squeeze(nanmean(w,3)),1)>0);
for i = active
   motif = squeeze(temp(:,i,:));
   motif = reshape(motif,[68,68,size(temp,3)]);
   for j = 1:size(temp,3)
      imagesc(motif(:,:,j),[0 0.75]); colorbar     
      axis equal; axis off
      title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
      pause();
   end      
end

%%
temp = NaN(numel(nanpxs)+size(W_basis,1),size(W_basis,2),size(W_basis,3)); 
goodpxl = 1:size(temp,1);
temp(~ismember(goodpxl,nanpxs),:,:)=W_basis;

figure; 
active = find(nanmean(squeeze(nanmean(W_basis,3)),1)>0);
for i = active
   motif = squeeze(temp(:,i,:));
   motif = reshape(motif,[68,38,size(temp,3)]);
   for j = 1:size(temp,3)
      imagesc(motif(:,:,j),[0 0.25]); colorbar     
      axis equal; axis off
      title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
      pause(0.01);
   end      
end


%% 4032 4046
figure;

% temp = conditionDffMat(data_norm',nanpxs,[],[68,38,size(data_norm,2)]);
% temp = data;
% temp = conditionDffMat(data',nanpxs,[],[68,38,size(data,2)]);
% temp = data_norm;
temp = data(:,:,10000:20000);
for j = 1:size(temp,3)
  imagesc(temp(:,:,j),[0 .25]);      
  axis equal; axis off
  colormap magma
  title(sprintf('motif %d, frame %d',i,j),'fontsize',12,'fontweight','normal');
  pause(0.01);
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