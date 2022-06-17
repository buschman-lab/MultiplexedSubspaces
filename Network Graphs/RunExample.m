%@Example
%Camden MacDowell - timeless

%% addpaths
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));

%% load fitted models (~60 secs)
data = LoadSubspaceData('in_grouped');

% figure parameters class to make life beautiful
fp = fig_params_cortdynamics;
fprintf('Data loaded.\n')

%% compile our interactions (~10 seconds) | USING A UNIFORM SAMPLE (same dimension)
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 4; %number of dimensions to use for exhaustive search
[rho, area_label] = CompileInteractions(data,useProj,ndim,'uniform');
fprintf('Interactions compiled.\n')

%% Plot example interactions 
%rho is rho{recording,motif}{area}(:,:,dimension)

%recording 1, motif 5; correlations between all pairs of target regions
%within the beta weights of area 2 (MOs), dimension 1
r = rho{1,5}{2}(:,:,1);
figure; hold on;
imagesc(r); c=colorbar; colormap parula; 
ylabel(c,'Similarity (rho)','FontSize',14);
set(gca,'xlim',[0.5,8.5],'ylim',[0.5,8.5],'xtick',[1:8],'ytick',[1:8],'xticklabel',area_label,'YTickLabel',area_label,'XTickLabelRotation',90)
%gut check is that MOs should be all NaNs here (since it can project to itself)
fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
fp.FormatAxes(gca);
title('Interactions within MOs dim 1','FontWeight','normal')

%recording 6, motif 2; correlations between all pairs of target regions
%within the beta weights of area 8 (VIS), dimension 3
r = rho{6,2}{8}(:,:,3);
figure; hold on;
imagesc(r); c=colorbar; colormap parula; 
ylabel(c,'Similarity (rho)','FontSize',14);
set(gca,'xlim',[0.5,8.5],'ylim',[0.5,8.5],'xtick',[1:8],'ytick',[1:8],'xticklabel',area_label,'YTickLabel',area_label,'XTickLabelRotation',90)
%gut check is that MOs should be all NaNs here (since it can project to itself)
fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
fp.FormatAxes(gca);
title('Interactions within VIS dim 3','FontWeight','normal')

%% Average across motifs and recordings (which is what we probably want to do for subsequent analysis)
rho_avg = cat(2,rho{:});
%fisher transform, average, fisher inverse
%note; you do get some structure difference when using mean versus median.
%Might be worth looking at both... 
rho_avg = arrayfun(@(n) fisherInverse(nanmedian(fisherZ(cat(4,rho_avg{n,:})),4)),1:size(rho_avg,1),'UniformOutput',0); 
%now you have an cell array containing the interactions within each region

%% Example agglomerative plots
cur_a = 1; %plot for HIPP
for cur_d = 1:4
   figure; hold on; 
   r = rho_avg{cur_a}(:,:,cur_d).^2; %plot the rsq for better line contrast in size
   r(isnan(r))=0;
   circularGraph(r,'colormap',repmat(fp.c_area(cur_a,:),8,1),'label',area_label,'normVal',0.25); %plot
   title(sprintf('%s Dim %d Subspace Networks (rsq range: %0.2f-%0.2f)',area_label{cur_a},cur_d,min(r(r>0)),max(r(:))),'FontWeight','normal')
   set(findall(gca, 'type', 'text'), 'visible', 'on') 
   fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
%saving function | need to supply a savedir
% saveCurFigs(get(groot,'Children'),{'-dpng'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png
% saveCurFigs(get(groot,'Children'),{'-dpng','-dsvg'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png and svg   
end

%% compile our interactions (~10 seconds) | USING AN Exhaustive serach (i.e, best correlation between two target areas within a source region but mix and matching dimensions
useProj = 1; %def = 1; use the timepoint x trial projection along betas. 0 = just correlate betas
ndim = 4; %number of dimensions to use for exhaustive search
[rho, area_label, pval] = CompileInteractions(data,useProj,ndim,'exhaustive');
fprintf('Exhaustive compilation of interactions complete.\n');

%% Plot example: 
%rho is now rho{recording,motif}{area}(ndim*narea,ndim*narea)
%recording 1, motif 5; correlations between all pairs of target regions
%within the beta weights of area 2 (MOs)
r = rho{1,5}{2}(:,:);
p = pval{1,5}{2}(:, :);

%plot all the correlations
figure; hold on;
imagesc(r); c=colorbar; colormap parula; 
ylabel(c,'Similarity (rho)','FontSize',14);
lab = repmat(area_label,1,ndim)'; 
lab = lab(:);
suffix = repmat(arrayfun(@(x) num2str(x), 1:ndim,'UniformOutput',0),numel(area_label),1)';
suffix = suffix(:); 
lab = cellfun(@(x,y) [x,' dim ', y],lab,suffix,'UniformOutput',0);
set(gca,'xlim',[0.5,size(r,1)+.5],'ylim',[0.5,size(r,1)+.5],'xtick',[1:size(r,1)],'ytick',[1:size(r,1)],'xticklabel',lab(:),'YTickLabel',lab,'XTickLabelRotation',90)
% fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
fp.FormatAxes(gca);
title('Interactions within MOs','FontWeight','normal')

%plot the p-value of the correlations
figure; hold on;
imagesc(-log10(p)); c=colorbar; colormap parula; 
ylabel(c,'-log10(p-Value)','FontSize',14);
lab = repmat(area_label,1,ndim)'; 
lab = lab(:);
suffix = repmat(arrayfun(@(x) num2str(x), 1:ndim,'UniformOutput',0),numel(area_label),1)';
suffix = suffix(:); 
lab = cellfun(@(x,y) [x,' dim ', y],lab,suffix,'UniformOutput',0);
set(gca,'xlim',[0.5,size(p,1)+.5],'ylim',[0.5,size(p,1)+.5],'xtick',[1:size(p,1)],'ytick',[1:size(p,1)],'xticklabel',lab(:),'YTickLabel',lab,'XTickLabelRotation',90)
% fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
fp.FormatAxes(gca);
title('Interactions within MOs','FontWeight','normal')

%% Plot the average
rho_avg = cat(2,rho{:});
%fisher transform, average, fisher inverse
%note; you do get some structure difference when using mean versus median.
%Might be worth looking at both... 
rho_avg = arrayfun(@(n) fisherInverse(nanmedian(fisherZ(cat(4,rho_avg{n,:})),4)),1:size(rho_avg,1),'UniformOutput',0); 
%now you have an cell array containing the interactions within each region

%% Example average plots (obviously you'll want to take some peak or average here ... the exhuastive has too much going on. (fyi you can change line thickness within the circular graph function lines 90-94
for cur_a = 1:8
   figure; hold on; 
   r = rho_avg{cur_a}.^2; %plot the rsq for better line contrast in size
   r(isnan(r))=0;
   circularGraph(r,'colormap',repmat(fp.c_area(cur_a,:),8*ndim,1),'label',lab,'normVal',0.5); %plot
   title(sprintf('%s Subspace Networks (rsq range: %0.2f-%0.2f)',area_label{cur_a},min(r(r>0)),max(r(:))),'FontWeight','normal')
   set(findall(gca, 'type', 'text'), 'visible', 'on') 
   fp.FigureSizing(gcf,[3 3 16 16],[5 5 24 24])
%saving function | need to supply a savedir
% saveCurFigs(get(groot,'Children'),{'-dpng'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png
% saveCurFigs(get(groot,'Children'),{'-dpng','-dsvg'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png and svg   
end

%% Thresholded connection plots

for cur_a = 1:8
   figure; hold on; 
   p = pval{1,1}{cur_a}; %plot the rsq for better line contrast in size
   p = (p <= 0.05/sum(sum(~isnan(p))));
   p(isnan(p))=1;
   circularGraph(p,'colormap',repmat(fp.c_area(cur_a,:),8*ndim,1),'label',lab,'normVal',1); %plot
   title(sprintf('%s Subspace Networks (rsq range: %0.2f-%0.2f)',area_label{cur_a},min(p(p>0)),max(p(:))),'FontWeight','normal')
   set(findall(gca, 'type', 'text'), 'visible', 'on') 
   fp.FigureSizing(gcf,[3 3 16 16],[5 5 24 24])
%saving function | need to supply a savedir
% saveCurFigs(get(groot,'Children'),{'-dpng'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png
% saveCurFigs(get(groot,'Children'),{'-dpng','-dsvg'},sprintf('SubspaceNetworks_%sdim%d',area_label{cur_a},cur_d),savedir,0); %for png and svg   
end

