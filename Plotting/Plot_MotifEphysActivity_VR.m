function Plot_MotifEphysActivity_VR(cur_rec,motif_num,fast_load,type,savedir)
if nargin <4; type = 'mean'; end
%Camden MacDowell - timeless

fp = fig_params_cortdynamics;
if fast_load==0
    [~,~,~,EphysPath,~] = LoadDataDirectories_VR(cur_rec);
    [st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
else %preprocessed
    [rec_name,~,~,EphysPath] = LoadDataDirectories_VR(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif1.mat'],0,{'\\cup\buschman\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA'}); %first motif has all that you need
    temp = cellfun(@(x) load(x,'st_mat','st_depth'),fn);  
    st_depth = temp.st_depth;
    st_mat = temp.st_mat; 
end
[~,ImgPath,~,~,motif_fits] = LoadDataDirectories_VR(cur_rec);
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct
%remove noise motif
motif_onset(fp.noisemotif) = [];

%load the imaging data
load(ImgPath,'data_norm','nanpxs');

neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

win = [-4,22]; %this matches what we use in the subspace analyses
%get the average response per area
for i = 1:numel(motif_num)
    [~,trig_st] = ParseByOnset([],st_mat,motif_onset,win,motif_num(i));
    %VR added format neu_area as array
    parse_neu_area = [neu_area{1,1},neu_area{1,2},neu_area{1,3},neu_area{1,4}];
    [area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general'); %replaced parse_neu_area w/ neu_area
    %remove the rare edge case where a motif begins at the start (no baseline)
    area_val = RemoveEdgeTrials(area_val);

    %clean up areas %third input is the min # of spikes to keep area
    [area_val, area_label] = CleanUpAreas(area_val, area_label, 10);
 
    %parse imaging by onset
    trig_dff = ParseByOnset(data_norm',[],motif_onset,win,motif_num(i));
    %get average activity across neurons and trials
%     ttvdff = squeeze(nanmean(trig_dff,1));
%     ttvdff = nanmean(abs(ttvdff-nanmean(ttvdff,2)),2);
%     ttvdff = (ttvdff-nanmin(ttvdff))/(max(ttvdff)-min(ttvdff));
    dff = squeeze(nanmean(trig_dff,[1,3]));
    dff = (dff-nanmin(dff))/(max(dff)-min(dff));
    

    switch type
        case 'mean'
            %loop through area
            x = win(1):win(2);
            figure; hold on; 
            plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)
            y_all = cell(1,numel(area_val));
            for j = 1:numel(area_val)           
                y = area_val{j}; %just plot this motif         
                %normalize to baseline
                y = normalizeToBaseline(y,[1:4],'mean');               
                %avg across trials and save off for raster 
                y_all{j} = y(:,5:end,:);
                y_all{j} = squeeze(nanmean(y_all{j},3));
                %average across neurons
                y = squeeze(nanmean(y,1)); 
                y = y(5:end,:);
                shadedErrorBar(x(5:end),nanmean(y,2)-mean(nanmean(y,2)),sem(y,2),'lineprops',{'color',[fp.c_area(j,:),0.25],'linewidth',2});
            end            
            fp.FormatAxes(gca);              
            xlim([0,x(end)])
            xvals = get(gca,'xtick');            
            xlabel('time relative to onset (ms)');
            ylabel('Firing Rate (Baseline Normalized)')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            legend(cat(1,{''},area_label),'location','bestoutside')
            fp.FigureSizing(gcf,[3 2 4 4],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            
            %organize activity within each area by peak activity            
            y = cellfun(@(yy) OrderPSTH(yy), y_all,'UniformOutput',0);
            y = cat(1,y{:});
            tempdff = (0.1*size(y,1)*dff)+size(y,1);%adjust to the size of the axis  
            tempdff = tempdff(5:end);
            %for visualization, mean center
            y = y-nanmean(y,2);
            cval = [0 prctile(y(:),99)];
            figure; hold on; imagesc(y,cval); colormap(flipud(gray));
            c=colorbar;
            AddAreaPSTH(area_val,area_label,0.5,0)
            plot(tempdff,'linewidth',2,'color','k')

            xlim([0.5,x(end)+1])
            ylim([0 size(y,1)+0.5])
            xvals = get(gca,'xtick');            
            xlabel('time relative to onset (ms)');
            ylabel(c,'Firing Rate (Baseline Normalized)')
            ylabel('neurons')
            set(gca,'YAxisLocation','right','Clipping','off')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            fp.FigureSizing(gcf,[5 2 3 6],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            fp.FormatAxes(gca)     
            

            %same but showing trial-to-trial variance
            x = win(1):win(2);
            figure; hold on; 
            plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)
            y_all = cell(1,numel(area_val));
            for j = 1:numel(area_val)           
                y = area_val{j}; %just plot this motif         
                %normalize to baseline
                y = normalizeToBaseline(y,[1:4],'mean');
                %get trail to trail variance
                y = abs(y-nanmean(y,3));
                %avg across trials and save off for raster 
                y_all{j} = y(:,5:end,:);
                y_all{j} = squeeze(nanmean(y_all{j},3));
                %average across neurons
                y = squeeze(nanmean(y,1));
                y = y(5:end,:);
                shadedErrorBar(x(5:end),nanmean(y,2)-mean(nanmean(y,2)),sem(y,2),'lineprops',{'color',[fp.c_area(j,:),0.25],'linewidth',2});
            end
            fp.FormatAxes(gca);  
            xlim([0,x(end)])
            xvals = get(gca,'xtick');
            xlabel('time relative to onset (ms)');
            ylabel('trial-to-trial |\sigma^2|')
            title(sprintf('Motif %d',i),'FontWeight','normal'); 
            legend(cat(1,{''},area_label),'location','bestoutside')
            fp.FigureSizing(gcf,[3 2 4 4],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            
            %organize activity within each area by peak activity            
            y = cellfun(@(yy) OrderPSTH(yy), y_all,'UniformOutput',0);
            y = cat(1,y{:});
            tempdff = (0.1*size(y,1)*dff)+size(y,1);%adjust to the size of the axis  
            tempdff = tempdff(5:end);
            %for visualization, mean center
            y = y-nanmean(y,2);
            cval = [0 prctile(y(:),99)];
            figure; hold on; imagesc(y,cval); colormap(flipud(gray));
            c=colorbar;
            AddAreaPSTH(area_val,area_label,0.5,0)
            plot(tempdff,'linewidth',2,'color','k')

            xlim([0.5,x(end)+1])
            ylim([0 size(y,1)+0.5])
            xvals = get(gca,'xtick');            
            xlabel('time relative to onset (ms)');
            ylabel(c,'trial-to-trial |\sigma^2| in FR')
            ylabel('neurons')
            set(gca,'YAxisLocation','right','Clipping','off')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            fp.FigureSizing(gcf,[5 2 3 6],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            fp.FormatAxes(gca) 



        case 'meansubtract'
            %loop through area
            x = win(1):win(2);
            figure; hold on; 
            plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)
            y_all = cell(1,numel(area_val));
            for j = 1:numel(area_val)           
                y = area_val{j}; %just plot this motif         
                %normalize to baseline
                y = normalizeToBaseline(y,[1:4],'meansubtract');
                %avg across trials and save off for raster 
                y_all{j} = y(:,1:end,:);
                y_all{j} = squeeze(nanmean(y_all{j},3));
                %average across neurons
                y = squeeze(nanmean(y,1)); 
                shadedErrorBar(x,nanmean(y,2),sem(y,2),'lineprops',{'color',[fp.c_area(j,:),0.25],'linewidth',2});
            end
            fp.FormatAxes(gca);  
            xlim([x(1),x(end)])
            xvals = get(gca,'xtick');
            xlabel('time relative to onset (ms)');
            ylabel('Firing Rate (Baseline Subtracted)')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            legend(cat(1,{''},area_label),'location','bestoutside')
            yval = get(gca,'ylim');
            plot([0,0],yval,'color','r','linestyle',':','linewidth',2)
            ylim(yval);     
            fp.FigureSizing(gcf,[3 2 4 4],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
    
            %organize activity within each area by peak activity            
            y = cellfun(@(yy) OrderPSTH(yy), y_all,'UniformOutput',0);
            y = cat(1,y{:});
            tempdff = (0.1*size(y,1)*dff)+size(y,1);%adjust to the size of the axis  
            %for visualization, mean center
            y = y-nanmean(y,2);
            cval = [0 prctile(y(:),99)];
            figure; hold on; imagesc(y,cval); colormap(flipud(gray));
            c=colorbar;
            AddAreaPSTH(area_val,area_label,0.5,0)
            plot(tempdff,'linewidth',2,'color','k')

            xlim([0.5,x(end)+1])
            ylim([0 size(y,1)+0.5])
            xvals = get(gca,'xtick');            
            xlabel('time relative to onset (ms)');
            ylabel(c,'trial-to-trial |\sigma^2| in FR')
            ylabel('neurons')
            set(gca,'YAxisLocation','right','Clipping','off')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            fp.FigureSizing(gcf,[5 2 3 6],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            fp.FormatAxes(gca) 

            %same but showing trial-to-trial variance
            x = win(1):win(2);
            figure; hold on; 
            plot([x(1),x(end)],[0 0],'color','k','linestyle',':','linewidth',2)
            y_all = cell(1,numel(area_val));
            for j = 1:numel(area_val)           
                y = area_val{j}; %just plot this motif         
                %normalize to baseline
                y = normalizeToBaseline(y,[1:4],'meansubtract');
                %get trail to trail variance
                y = abs(y-nanmean(y,3));
                %avg across trials and save off for raster 
                y_all{j} = y(:,1:end,:);
                y_all{j} = squeeze(nanmean(y_all{j},3));
                %average across neurons
                y = squeeze(nanmean(y,1));
                y(1:4,:)=NaN;
                shadedErrorBar(x,nanmean(y,2)-min(nanmean(y,2)),sem(y,2),'lineprops',{'color',[fp.c_area(j,:),0.25],'linewidth',2});
            end
            fp.FormatAxes(gca);  
            xlim([x(1),x(end)])
            xvals = get(gca,'xtick');
            xlabel('time relative to onset (ms)');
            ylabel('trial-to-trial |\sigma^2|')
            title(sprintf('Motif %d',i),'FontWeight','normal'); 
            legend(cat(1,{''},area_label),'location','bestoutside')
            yval = get(gca,'ylim');
            plot([0,0],yval,'color','r','linestyle',':','linewidth',2)
            ylim(yval); 
            fp.FigureSizing(gcf,[3 2 4 4],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));

            %organize activity within each area by peak activity            
            y = cellfun(@(yy) OrderPSTH(yy), y_all,'UniformOutput',0);
            y = cat(1,y{:});
            tempdff = (0.1*size(y,1)*dff)+size(y,1);%adjust to the size of the axis  
            %for visualization, mean center
            y = y-nanmean(y,2);
            cval = [0 prctile(y(:),99)];
            figure; hold on; imagesc(y,cval); colormap(flipud(gray));
            c=colorbar;
            AddAreaPSTH(area_val,area_label,0.5,0)
            plot(tempdff,'linewidth',2,'color','k')

            xlim([0.5,x(end)+1])
            ylim([0 size(y,1)+0.5])
            xvals = get(gca,'xtick');            
            xlabel('time relative to onset (ms)');
            ylabel(c,'trial-to-trial |\sigma^2| in FR')
            ylabel('neurons')
            set(gca,'YAxisLocation','right','Clipping','off')
            title(sprintf('Motif %d',motif_num(i)),'FontWeight','normal'); 
            fp.FigureSizing(gcf,[5 2 3 6],[2 10 12 12])
            set(gca,'xtick',xvals,'XTickLabel',round(xvals*(1000/15)));
            fp.FormatAxes(gca) 

    end
    set(gca,'clim',[0 0.1])
    saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},sprintf('Motif_triggeredActivity_motif%d_rec%d',i,cur_rec),savedir,0); close all
end



end










