function WidefieldProbeInsertionLocations(recording_paths,labels,save_dir)
%Camden MacDowell - timeless

if nargin <1
   recording_paths = {'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_11_2021_1dff_combined.mat',...
        'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse331_RestingState_NP_06_12_2021_1dff_combined.mat',...
        'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_07_2021_1dff_combined.mat',...
        'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse332_RestingState_NP_06_08_2021_1dff_combined.mat',...
        'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_09_2021_1dff_combined.mat',...
        'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\Mouse334_RestingState_NP_06_10_2021_1dff_combined.mat'};
end

if nargin <2
   labels = {'331-1','331-2','332-1','332-2','334-1','334-2'}';
end

if nargin <3
    save_dir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ProbeLocalization';
end

%% meat and potatoes
close all;
%Bregma coordinates
BregmaX = 1.97; %conversion ratio for finding location to plot bregma
BregmaY = 2.27; %conversion ratio for finding locatio to plot bremga
bX = 68/BregmaX+1;
bY = 68/BregmaY;

%load the probe coordinates from the widefield imaging
p = cellfun(@(x) load(x,'probe_coords','opts'), recording_paths,'UniformOutput',1);
opts = p(1).opts;

all_mask = arrayfun(@(n) p(n).opts,1:6,'UniformOutput',1);
all_mask = sum(cat(3,all_mask(:).mask),3)==6;

%plot the coordinates of the insertions
col = distinguishable_colors(6);
%colored same per mouse
col(2:2:end,:) = col(1:2:end,:);
markertype = repmat({'x','*'},1,3);
figure; hold on;
imagesc(imresize(opts.mask,1/opts.spatial_bin_factor))
colormap gray
axis off
pp = cell(1,numel(recording_paths));
for i = 1:numel(recording_paths)
    pp{i} = arrayfun(@(n) plot(p(i).probe_coords{n}(1,1),p(i).probe_coords{n}(1,2),'linewidth',1,'color',col(i,:),'marker',markertype{i},'markersize',10,'linestyle','none'),1:numel(p(1).probe_coords),'UniformOutput',0);
end
set(gca,'ydir','reverse')
pp = cellfun(@(x) x{1},pp,'UniformOutput',1);
legend(pp,labels,'Location','bestoutside')

%Addbregma
scatter(bX,bY,600,'.','MarkerFaceColor',[1 0.2 .2],'MarkerEdgeColor',[1 0.2 .2]); 
ylim([0 68]);
xlim([0 68]);

%Set Axes
ylim([0 68]);
xlim([0 68]);
axis square

%Using the combined mask across animals
col = distinguishable_colors(6);
%colored same per mouse
col(2:2:end,:) = col(1:2:end,:);
markertype = repmat({'x','*'},1,3);
figure; hold on;
imagesc(imresize(all_mask,1/opts.spatial_bin_factor))
colormap gray
axis off
pp = cell(1,numel(recording_paths));
for i = 1:numel(recording_paths)
    pp{i} = arrayfun(@(n) plot(p(i).probe_coords{n}(1,1),p(i).probe_coords{n}(1,2),'linewidth',1,'color',col(i,:),'marker',markertype{i},'markersize',10,'linestyle','none'),1:numel(p(1).probe_coords),'UniformOutput',0);        
end
set(gca,'ydir','reverse')
title('combined mask');
pp = cellfun(@(x) x{1},pp,'UniformOutput',1);
legend(pp,labels,'Location','bestoutside')

%Addbregma
scatter(bX,bY,600,'.','MarkerFaceColor',[1 0.2 .2],'MarkerEdgeColor',[1 0.2 .2]); 
ylim([0 68]);
xlim([0 68]);

%Set Axes
ylim([0 68]);
xlim([0 68]);
axis square


%save off
saveCurFigs(get(groot, 'Children'),{'-dsvg','-dpng'},'ProbeInsertionLocations',save_dir,0); close all


end %function













