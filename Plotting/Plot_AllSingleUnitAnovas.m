function Plot_AllSingleUnitAnovas(savedir)
%Camden - timeless
if nargin<1; savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SingleUnitAnovas'; end
if ~exist(savedir,'dir'); mkdir(savedir); end
    
%Loop through each recording
for cur_rec = 1:6
    fprintf('\n\tworking on rec %d',cur_rec);
    win = [-5 15];
    [rec_name,~,~,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);
    [sig_motif,weight_motif,pref_motif,neu_area,discrim_mat,num_discrim] = SingleUnitANOVA(EphysPath,motif_fits,win);

    
    %remove the tiny number of neurons in white matter or unlabeled regions
    bad_idx = cellfun(@(x) find(strcmp([neu_area.parent_label],x)==1),{'cc','fiber tracts','na'},'UniformOutput',0);
    bad_idx = [bad_idx{:}];
    
    neu_area(bad_idx) = [];
    sig_motif(bad_idx,:) = [];
    weight_motif(bad_idx,:) = [];
    pref_motif(bad_idx,:) = [];
    discrim_mat(bad_idx,:,:) = [];
    discrim_mat(:,bad_idx,:) = [];
    num_discrim(bad_idx,:)=[];

    figure('position',[684    68   339   928]); hold on;
    %anything below our significance threshold (per timepoint) denote as NaN
    sig_map = sig_motif; 
    %threshold by all neurons and all timepoints (probably not the write way)
%     sig_map(sig_map>=0.05/numel(sig_map))=NaN;
%     sig_map(sig_map<0.05/numel(sig_map))=1;
    %threshold significance per neuron
    sig_map(sig_map>=0.05/size(sig_map,2))=NaN;
    sig_map(sig_map<0.05/size(sig_map,2))=1;
    
    
    %colormap per motif
    motif_list = [1,3,4,5,6,7,8,9,10,11,12,13,14,15];
    col = getColorPalet(numel(motif_list));
    col = cat(1,[1 1 1],col);
    temp = pref_motif;
    temp(isnan(sig_map))=0;
    imagesc(temp); colormap(col); hold on; 
    ylim([0.5,size(sig_map,1)]);
    set(gca,'ydir','reverse')
    c = getColorPalet(30);
    PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
    set(gca,'ytick','')
    line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset   

    %loop through each motif    
    for i = 1:numel(motif_list)
        figure('position',[684    68   339   928]); hold on; 
        temp = pref_motif;
        temp(isnan(sig_map))=0; 
        temp(temp~=i)=0;
        imagesc(temp); colormap(cat(1,col(1,:),col(i,:))); hold on; 
        ylim([0.5,size(sig_map,1)]);
        set(gca,'ydir','reverse')
        c = getColorPalet(30);
        PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
        set(gca,'ytick','')
        line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset  
        title(sprintf('preference for motif %d',motif_list(i)),'fontweight','normal')
    end
    saveCurFigs(get(groot, 'Children'),{'-dpng'},rec_name,savedir,0); close all 
    
    %Also plot a bar chart showing the relative mixing within the region
    %loop through all areas    
    area_list = [neu_area.parent_label];
    unique_area = unique(area_list);
    counts = NaN(numel(unique_area),numel(motif_list));
    for i = 1:numel(unique_area)
        idx = strcmp(area_list,unique_area(i));         
        temp = pref_motif(idx,:);
        temp(isnan(sig_map(idx,:)))=0; 
        counts(i,:) = arrayfun(@(n) sum(temp(:)==n), motif_list,'UniformOutput',1);                
    end
    
    subcort_list = {'EPI','ACAd','ILA','MED','DG','CA','PL','fxs','ifbst','mfbse'};
    figure; hold on; 
    [num_rows, num_col]=numSubplot(size(counts,1),0.5);
    for i = 1:size(counts,1)       
       subplot(num_rows,num_col,i); hold on; 
       bar(counts(i,:),'EdgeAlpha',0,'FaceColor',[0.5 0.2 0.2],'FaceAlpha',0.5)
       if ismember(unique_area{i},subcort_list)
           title(unique_area{i},'FontWeight','bold')
       else
           title(unique_area{i},'FontWeight','normal')
       end
       set(gca,'xtick',1:numel(motif_list),'XTickLabel',motif_list)
    end
        
    saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name,'mixingwithinarea'],savedir,0); close all 
    
    %get the average within the brain areas 
    within_ratio = NaN(size(discrim_mat,1),size(discrim_mat,3));
    for cur_t = 1:size(discrim_mat,3)
        for cur_n = 1:size(discrim_mat,1)
            temp = squeeze(discrim_mat(cur_n,:,cur_t))'; 
            %set self to NaN
            temp(cur_n,1)=NaN;
            %split by area
            [area_val, area_lbl] = ParseByArea(temp,neu_area,'parent');
            %get within similarity (could also normalize by across)
            idx = strcmp(area_lbl,neu_area(cur_n).parent_label)==1;
            within_ratio(cur_n,cur_t) = nanmean(area_val{idx});        
        end
    end        
            
    %Plot the similarity in discriminatability by region 
    figure('position',[684    68   500   928]); hold on; 
    imagesc(within_ratio,[0 0.5]); colormap(flipud(gray)); hold on; 
    ylim([0.5,size(within_ratio,1)]);
    set(gca,'ydir','reverse')
    c = getColorPalet(30);
    PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
    set(gca,'ytick','')
    line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset       
    colorbar;
    title({'Ratio of shared pairs of','discriminatble motifs with','units from same area'});
    saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name,'discrim_shared'],savedir,0); close all 

    figure('position',[684    68   500   928]); hold on; 
    imagesc(num_discrim,[0 20]); colormap(flipud(gray)); hold on; 
    ylim([0.5,size(num_discrim,1)]);
    set(gca,'ydir','reverse')
    c = getColorPalet(30);
    PlotProbeAnatomy(gca, neu_area, 0, 'parent',1,0,c(randperm(size(c,1),size(c,1)),:));
    set(gca,'ytick','')
    line(repmat(abs(win(1))+1,2),get(gca,'ylim'),'color','r','linestyle',':','linewidth',2); %line showing onset       
    colorbar;
    title({'Number of Disciminitable','Pairs of Motifs'});
    
    saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name,'discrim_number'],savedir,0); close all 
    
    %average per brain ares
    [area_val, area_lbl] = ParseByArea(within_ratio,neu_area,'parent');
    area_val = cellfun(@(x) nanmean(x(:)),area_val,'UniformOutput',1);
    figure; hold on; 
    bar(area_val*100)
    set(gca,'xtick',1:size(area_lbl),'XTickLabel',area_lbl,'XTickLabelRotation',45);
    ylabel('percentage of shared discimination')
    
    %average per brain ares
    [area_val, area_lbl] = ParseByArea(num_discrim,neu_area,'parent');
    area_val = cellfun(@(x) nanmean(x(:)),area_val,'UniformOutput',1);
    figure; hold on; 
    bar(area_val)
    set(gca,'xtick',1:size(area_lbl),'XTickLabel',area_lbl,'XTickLabelRotation',45);
    ylabel('number of pairs of motifs seperable') 
    
        
    saveCurFigs(get(groot, 'Children'),{'-dpng'},[rec_name,'discrim_by_area'],savedir,0); close all 
 
    
end %rec loop

end %function end















