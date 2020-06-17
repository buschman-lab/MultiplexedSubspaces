function ImpactOfHemoCorrection(dff,dff_b,dff_h)
%Camden MacDowell - timeless. Stack is a x,y,time imaging dff (e.g. post
%hemo correction). This looks at the impact on pca denoising

if size(dff,3)>7500 %trim is super big
    dff = dff(:,:,end-7500:end);
    dff_b = dff_b(:,:,end-7500:end);
    dff_h = dff_h(:,:,end-7500:end);
end

[x,y,~] = size(dff);

%reg dff
fprintf('\n Working on dff');
[dff_flat,nanpxs] = conditionDffMat(dff); 
[~, score, ~, ~, ~, ~] = pca(dff_flat');
score = conditionDffMat(score',nanpxs,[],[x,y,size(score,2)]);

%make figure;
num_pcs = 10;
figure('position',[223,558,1546,420]); hold on;
for i = 1:num_pcs  
    subplot(1,num_pcs,i);
    set(gca,'ydir','reverse');
    imagesc(score(:,:,i));
    axis equal
end
sgtitle('Corrected PCs');


%blue dff
fprintf('\n Working on dff_b');
[dff_flat,nanpxs] = conditionDffMat(dff_b); 
[~, score, ~, ~, ~, ~] = pca(dff_flat');
score = conditionDffMat(score',nanpxs,[],[x,y,size(score,2)]);

%make figure;
num_pcs = 10;
figure('position',[223,558,1546,420]); hold on;
for i = 1:num_pcs  
    subplot(1,num_pcs,i);
    set(gca,'ydir','reverse');
    imagesc(score(:,:,i));
    axis equal
end
sgtitle('Dff_b PCs');





%hemo dff
fprintf('\n Working on dff_h');
[dff_flat,nanpxs] = conditionDffMat(dff_h); 
[~, score, ~, ~, ~, ~] = pca(dff_flat');
score = conditionDffMat(score',nanpxs,[],[x,y,size(score,2)]);

%make figure;
num_pcs = 10;
figure('position',[223,558,1546,420]); hold on;
for i = 1:num_pcs  
    subplot(1,num_pcs,i);
    set(gca,'ydir','reverse');
    imagesc(score(:,:,i));
    axis equal
end
sgtitle('Dff_h PCs');
end
