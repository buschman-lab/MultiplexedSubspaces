function handles = AnalyzeResiduals(X,Xhat,nanpxs,save_dir)
%Camden MacDowell. Looks for structure in residuals of CNMF fit
%X and Xhat are pxl x time matrices

if nargin <3; nanpxs = []; end
if nargin <4; save_dir = []; end

gp = general_params;

residuals = X-Xhat;

%Compare all 3
figure('position',[158 558 1645 420]); 
subplot(1,3,1); imagesc(X); title('original','fontweight','normal'); colorbar
subplot(1,3,2); imagesc(Xhat); title('reconstruction','fontweight','normal'); colorbar
subplot(1,3,3); imagesc(residuals); title('residuals','fontweight','normal'); colorbar
cmap = flipud(gray);
colormap(cmap)

%perform PCA to analyze spatial structure of residuals
[~, score_resid, ~, ~, explained, ~] = pca(residuals);
[~, score_x] = pca(X);

%plot variance explained
figure(); plot(cumsum(explained)); ylim([0, 100]);

num_comp = 20; %nuber of components to plot, 20 is fine
if ~isempty(nanpxs)%recondition matrices     
    score_resid = score_resid(:,1:num_comp);
    score_x = score_x(:,1:num_comp);
    score_resid = conditionDffMat(score_resid',nanpxs,[],[gp.pixel_dim,num_comp]);    
    score_x = conditionDffMat(score_x',nanpxs,[],[gp.pixel_dim,num_comp]);    
end

%plot the first 20 components of residuals
[nR,nC] = numSubplot(20,.5);
figure('units','normalized','position',[0 0 1 1]); hold on; 
for i = 1:num_comp
   subplot(nR,nC,i);
   imagesc(score_resid(:,:,i));
   axis equal; axis off;
end
title('Top 20 Spatial Components of Residuals');

%plot the first 20 components of original
[nR,nC] = numSubplot(20,.5);
figure('units','normalized','position',[0 0 1 1]); hold on; 
for i = 1:num_comp
   subplot(nR,nC,i);
   imagesc(score_x(:,:,i));
   axis equal; axis off;
end
title('Spatial Components of Original Data');

handles = get(groot, 'Children');

if ~isempty(save_dir)
    saveCurFigs(handles,'-dpng',sprintf('ResidualAnalysis'),save_dir,0); 
end

end



    