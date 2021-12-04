function RecordingAllignment()
%Camden MacDowell - timeless
%Figures to evaluate alignment between recordings

fp = fig_params;
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RecordingAllignment';
[dff_list,~] = GrabFiles('\w*_processed.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\DeconvolvedImaging'}); %select the preprocessed data (not the '_processed');
dff = cellfun(@(x) load(x,'data_norm','nanpxs'),dff_list,'UniformOutput',0);

%% Compare the allignment using 4 FC modules | Use non-craniotomy side
g = [18, 23; 37, 26; 41, 12; 51, 20; 48, 30];
%Get the seed corr map for each location for each rec
rho_map = NaN(68,68,size(g,1),numel(dff_list));
for cur_r = 1:numel(dff_list)
    data = conditionDffMat(dff{cur_r}.data_norm',dff{cur_r}.nanpxs);
    for i = 1:size(g,1)
       seed_vec = squeeze(nanmean(data(g(i,1)-1:g(i,1)+1,g(i,2)-1:g(i,2)+1,:),[1,2]));
       rho_map(:,:,i,cur_r) = conditionDffMat(corr(dff{cur_r}.data_norm'-nanmean(dff{cur_r}.data_norm'),seed_vec-nanmean(seed_vec))',dff{cur_r}.nanpxs);
    end
end

%keep no cranio hemi
rho_map_hemi = rho_map;
rho_map_hemi(:,34:end,:,:)=NaN;

%% Figure showing the alligment for each recording and each mouse
close all
for i = 1:size(rho_map,4) %animal loop
    [~,fn] = fileparts(dff_list{i});
    for j = 1:size(rho_map,3) %fc maps loop
        figure; hold on; 
        temp = squeeze(rho_map_hemi(:,:,j,i));
        temp(isnan(temp))=0;
        PlotMesoFrame(temp)
        caxis([0 0.75]); colormap gray; colorbar   
        scatter(g(j,2),g(j,1),80,'MarkerEdgeColor','m','LineWidth',2)
        set(gca,'ydir','reverse','xdir','reverse')        
        title(sprintf('FC %d',j))
    end
    saveCurFigs(get(groot, 'Children'),'-dpng',[fn,'_fcmap'],savedir,0); close all
end

%% Get the spatial correlation across animals
figure; hold on; 
% within_coords = {[1,2],[3,4],[5,6]}; %comparisons within animals
% within = cell(1,size(rho_map,3));
between_coords = {[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6],[3,5],[3,6],[4,5],[4,6]}; %comparisons across animals
between = cell(1,size(rho_map,3));
for i = 1:size(rho_map,3) %loop through each fc map
    temp = squeeze(rho_map_hemi(:,:,i,:));
    temp = conditionDffMat(temp);
    rho = corr(temp',temp','rows','complete');
    between{i} = cellfun(@(x) rho(x(1),x(2)),between_coords,'UniformOutput',1);
end
%violin plot comparing the location
CompareViolins([between{:}],fp);
ylim([0.5 1])
set(gca,'ytick',[0.5 0.75 1],'xtick',[])
ylabel('Spatial correlation')
fp.SetTitle(gca,{'Similarity',' across mice'});
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 2 4],[])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SpatialCorrelation',savedir,0); close all

%% Load retinotopy, manually register to the cranio | binarize and plot
%tbd



%% Retinotopy spatial correlation across animals vs within animals 




end