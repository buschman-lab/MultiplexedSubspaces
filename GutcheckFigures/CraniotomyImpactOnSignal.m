function CraniotomyImpactOnSignal(stack_path,processed_path)
%Camden MacDowell - timeless
%stack_path directs to a preprocessed (registered, spatially binned, stack)
%stack is NOT dff. 
%processed_path directs to the combined file (fully processed) with probe
%coordinates

%% compile data
if nargin <1
    stack_path = ['Z:\Rodent Data\Wide Field Microscopy\Neuropixels_Widefield_CorticalDynamics\',...
        'RestingState_Neuropixels\Mouse332_06_07_2021\Mouse332_RestingState_NP_06_07_2021_1\',...
        'Mouse332_RestingState_NP_06_07_2021_1_MMStack_Pos0_1.ome_stack.mat'];
end

if nargin <2
   processed_path = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined.mat';
end

stack = load(stack_path);
opts = stack.opts;
stack = stack.stack;
stack = filterstack(stack, 30, gp.w_filter_freq, gp.w_filter_type, 1, 0);
probe_coords = load(processed_path,'probe_coords');
probe_coords = probe_coords.probe_coords;

mask = repmat(imresize(opts.mask,[68 68]),1,1,size(stack,3));
stack(mask==0)=0;

opts.method = 'movingavg'; opts.method_window = 120;
dff = makeDFF(stack, opts, 'dff', opts.method_window);
dff(isnan(dff))=0;
for i = 1:size(dff,3); dff(:,:,i) = imgaussfilt(dff(:,:,i),'filtersize',[5 5]); end
dff(mask==0)=NaN;
opts.method = 'zscore'; opts.method_window = 120;
dfs = makeDFF(stack, opts, 'dff', opts.method_window);
dfs(isnan(dfs))=0;
for i = 1:size(dfs,3); dfs(:,:,i) = imgaussfilt(dfs(:,:,i),'filtersize',[5 5]); end
dfs(mask==0)=NaN;
stack(mask==0)=NaN;

%% compare raw, dff, dfs, in spatial view
close all; figure('position',[681 145 1065 834]); 
subplot(3,3,1); set(gca,'ydir','reverse'); hold on;  imagesc(nanmean(stack,3)); axis off; axis square; colorbar; title('mean')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,2); set(gca,'ydir','reverse'); hold on;  imagesc(nanmax(stack,[],3)); axis off; axis square; colorbar; title('max')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,3); set(gca,'ydir','reverse'); hold on;  imagesc(nanstd(stack,[],3)); axis off; axis square; colorbar; title('std')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)

%% dff
subplot(3,3,4); set(gca,'ydir','reverse'); hold on;  imagesc(nanmean(dff,3)); axis off; axis square; colorbar; title('mean')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,5); set(gca,'ydir','reverse'); hold on;  imagesc(nanmax(dff,[],3)); axis off; axis square; colorbar; title('max')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,6); set(gca,'ydir','reverse'); hold on;  imagesc(nanstd(dff,[],3)); axis off; axis square; colorbar; title('std')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)

%% dfs
subplot(3,3,7); set(gca,'ydir','reverse'); hold on;  imagesc(nanmean(dfs,3)); axis off; axis square; colorbar; title('mean')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,8); set(gca,'ydir','reverse'); hold on;  imagesc(nanmax(dfs,[],3)); axis off; axis square; colorbar; title('max')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)
subplot(3,3,9); set(gca,'ydir','reverse'); hold on;  imagesc(nanstd(dfs,[],3)); axis off; axis square; colorbar; title('std')
% cellfun(@(x) plot(x(1,1),x(1,2),'marker','.','color','r','markersize',10,'linewidth',2),probe_coords,'UniformOutput',0)


%% compare the distribution of pixel values 
%one hemi versus the other
L_mask = mask(:,:,1); L_mask(:,34:end) = 0;  
R_mask = mask(:,:,1); R_mask(:,1:34) = 0; %craniotomy side (not flipped)
c_rad = [-5:-1,1:5]; %center radius
s_rad = [-14:-8,8:14]; %surround radius
v_center = cellfun(@(x) cat(1,x(1,1)+c_rad, x(1,2)+c_rad),probe_coords,'UniformOutput',0); %rows are x and y
v_surround = cellfun(@(x) cat(1,x(1,1)+s_rad, x(1,2)+s_rad),probe_coords,'UniformOutput',0); %rows are x and y

%compare the methods
figure; 
subplot(3,1,1);hold on; 
raw_cent = stack(repmat(R_mask,1,1,size(stack,3))==1);
raw_surround = stack(repmat(L_mask,1,1,size(stack,3))==1);
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])
 
subplot(3,1,2);hold on; 
raw_cent = dff(repmat(R_mask,1,1,size(stack,3))==1);
raw_surround = dff(repmat(L_mask,1,1,size(stack,3))==1);
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])


subplot(3,1,3);hold on; 
raw_cent = dfs(repmat(R_mask,1,1,size(stack,3))==1);
raw_surround = dfs(repmat(L_mask,1,1,size(stack,3))==1);
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])

% center surround
figure; 
subplot(3,1,1); hold on; 
raw_cent = cellfun(@(x) stack(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
raw_cent = cat(3,raw_cent{:}); %combine 
raw_surround = cellfun(@(x) stack(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
raw_surround = cat(3,raw_surround{:}); %combine 
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])

 
subplot(3,1,2);hold on; 
raw_cent = cellfun(@(x) dff(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
raw_cent = cat(3,raw_cent{:}); %combine 
raw_surround = cellfun(@(x) dff(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
raw_surround = cat(3,raw_surround{:}); %combine 
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])

subplot(3,1,3);hold on; 
raw_cent = cellfun(@(x) dfs(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
raw_cent = cat(3,raw_cent{:}); %combine 
raw_surround = cellfun(@(x) dfs(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
raw_surround = cat(3,raw_surround{:}); %combine 
[s_counts,edges] = histcounts(raw_surround(:),100);s_counts = s_counts/sum(s_counts);
[c_counts,~] = histcounts(raw_cent(:),'BinEdges',edges);c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])






end


% 
% % center surround
% figure; 
% subplot(3,1,1); hold on; 
% raw_cent = cellfun(@(x) stack(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
% raw_cent = cat(3,raw_cent{:}); %combine 
% raw_surround = cellfun(@(x) stack(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
% raw_surround = cat(3,raw_surround{:}); %combine 
% [counts,edges] = histcounts(raw_surround(:),100);
% histogram(raw_cent(:),'BinEdges',edges,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
% histogram(raw_surround(:),'BinEdges',edges,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.75,'EdgeAlpha',0); 
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.5 0.5 0.5])
% [counts,~] = histcounts(raw_cent(:),'BinEdges',edges);
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
% 
%  
% subplot(3,1,2);hold on; 
% raw_cent = cellfun(@(x) dff(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
% raw_cent = cat(3,raw_cent{:}); %combine 
% raw_surround = cellfun(@(x) dff(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
% raw_surround = cat(3,raw_surround{:}); %combine 
% [counts,edges] = histcounts(raw_surround(:),100);
% histogram(raw_cent(:),'BinEdges',edges,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
% histogram(raw_surround(:),'BinEdges',edges,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.75,'EdgeAlpha',0); 
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.5 0.5 0.5])
% [counts,~] = histcounts(raw_cent(:),'BinEdges',edges);
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
% 
% subplot(3,1,3);hold on; 
% raw_cent = cellfun(@(x) dfs(x(1,:),x(2,:),:),v_center,'UniformOutput',0);
% raw_cent = cat(3,raw_cent{:}); %combine 
% raw_surround = cellfun(@(x) dfs(x(1,:),x(2,:),:),v_surround,'UniformOutput',0);
% raw_surround = cat(3,raw_surround{:}); %combine 
% [counts,edges] = histcounts(raw_surround(:),100);
% histogram(raw_cent(:),'BinEdges',edges,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
% histogram(raw_surround(:),'BinEdges',edges,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.75,'EdgeAlpha',0); 
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.5 0.5 0.5])
% [counts,~] = histcounts(raw_cent(:),'BinEdges',edges);
% plot(edges(1:end-1)+diff(edges)/2,smooth(counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])














