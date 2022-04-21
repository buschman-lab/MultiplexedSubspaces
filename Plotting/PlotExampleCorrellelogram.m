function PlotExampleCorrellelogram(cur_rec)
fp = fig_params_cortdynamics;
%Camden MacDowell - timeless
[area_val,area_label] = GetExampleSpikingData(cur_rec,1);

x = cat(1,area_val{:});
x = corr(x');
col = fp.c_area;

%%
figure; hold on;
imagesc(x,[0 0.04]); colorbar; 
colormap 'magma'
xlim([-5,size(x,2)+0.5])
ylim([0,size(x,1)+0.5])

for i = 1:numel(area_label)    
    if i>1
        temp = cat(1,area_val{1:i-1});
        idx = [size(temp,1)+1,size(temp,1)+size(area_val{i},1)];
    else
        idx = [1,size(area_val{i},1)];
    end
    plot([-5 -5],idx,'color',col(i,:),'linewidth',3); 
    text(-10, idx(1)+(idx(2)-idx(1))/2, area_label{i},'FontWeight','bold','HorizontalAlignment','right','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    plot(idx,[-5 -5],'color',col(i,:),'linewidth',3); 
    text(idx(1)+(idx(2)-idx(1))/2, -10, area_label{i},'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end
xlim([-15,size(x,2)+0.5])
ylim([-15,size(x,1)+0.5])
set(gca,'ytick',[],'xtick',[]);
fp.FormatAxes(gca); box on
set(gca,'Clipping','off')
axis square
end

