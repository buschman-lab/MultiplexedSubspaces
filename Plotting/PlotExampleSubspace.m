function PlotExampleSubspace(data,area_name,motif,cur_rec)
fp = fig_params_cortdynamics;

ndim = 15;
thresh = 0.80;

[rrr_mdl,area_label] = LoadVariable(data,'rel_performance',area_name);

%plot the relative activity 
col = fp.c_area; col = col(strcmp(area_label,area_name),:);
figure; hold on;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
if numel(motif)==1
    rrr_mdl = squeeze(rrr_mdl(:,motif,1:ndim));
    arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',1.5,'markersize',4), 1:size(rrr_mdl,1));
else
%     rrr_mdl = squeeze(nanmean(rrr_mdl,:)); %average across recordings
    rrr_mdl = reshape(rrr_mdl(:,motif,1:ndim),size(rrr_mdl,1)*numel(motif),ndim);
    arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.5,'markersize',4), 1:size(rrr_mdl,1));
end

plot(1:size(rrr_mdl,2), nanmean(rrr_mdl),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');  
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel({'Performance','(r^2 of explainable variance)'});
title(sprintf('%s | motif %d',area_name,motif),'fontweight','normal')
fp.FormatAxes(gca);  box on;
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])

%skip the rest if 

if numel(motif)==1
    %plot the current induction normalized per column
    ndim = 5; 

    b = cell(1,ndim);
    for i = 1:ndim
        b{i} = LoadVariable(data,'rrr_synapticweight_sorted',area_name,i);
    end
    %plot the average for that motif
    allbetas = cat(4,b{:});
    allbetas = squeeze(allbetas(cur_rec,motif,:,:));
    for i = 1:ndim
        temp = squeeze(allbetas(:,i));
        if sum(temp>0)<sum(temp<0)
            temp = temp*-1;
        end
        allbetas(:,i) = temp;
    end

    %uncomment to normalize by column
    % allbetas = allbetas./max(allbetas);
    allbetas(sum(isnan(allbetas),2)==ndim,:)=[];

    figure; hold on; 
    cval = prctile(allbetas(:),95);
    imagesc(allbetas,[-1*cval cval]); colormap(flipud(redgreencmap))
    xlim([0.5 ndim+0.5]);
    ylim([0.5 size(allbetas,1)]);
    c=colorbar;
    ylabel(c,'Beta weighted by FR var');

    reg = LoadVariable(data,'beta_region',area_name);
    lbl = squeeze(reg(1,motif,:));
    lblu = unique(lbl(~isnan(lbl)));
    c = getColorPalet(numel(lblu));
    p = cell(1,numel(lblu));
    xval = NaN(1,numel(lblu));
    for i = 1:numel(lblu)
       p{i} = plot(0.5*ones(1,sum(lbl==lblu(i))),find(lbl==lblu(i)),'color',c(i,:),'linewidth',2);
       xval(i) = mean(find(lbl==lblu(i)));
    end
    ylabel('brain area')
    xlabel('subspace dimension')
    xlim([0.5 ndim+0.5]);
    set(gca,'ytick',xval,'YTickLabel',area_label,'xtick',1:2:ndim);
    fp.FormatAxes(gca);  box on;
    fp.FigureSizing(gcf,[3 2 2 4],[10 10 20 10])
end
end %function end




















