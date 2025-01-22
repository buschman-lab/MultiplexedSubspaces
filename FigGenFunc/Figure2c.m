function Figure2c(b,area,area_label)
fp = fig_params_cortdynamics;
figure; hold on; 
bar(b,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'BarWidth',1,'FaceAlpha',1)
plot(1:numel(b),b,'linewidth',1.5,'color','r')
plot([-20,numel(b)+20],[0,0],'linewidth',1.5,'color','k')
fp.FormatAxes(gca); box on; grid on; 
set(gca,'ydir','reverse','YAxisLocation','right')
xlabel('neurons');
ylabel('subspace betas');
title({'visual subspace','dimension 1'},'fontweight','normal');
xlim([1,numel(b)])
camroll(-90)
fp.FigureSizing(gcf,[3 2 3 4],[2 10 15 15])

%plot the identity on top
id = unique(area);
for i = 1:7
    figure; hold on; 
    [y,edges] = histcounts(find(area==id(i)),0:20:numel(area)+20);   
    %interpolate
    y = interp1(edges(2:end),y,1:numel(area));
    y = smoothdata(y,'gaussian',50);
    %append the last data point for plotting
    y(isnan(y))=0;
    y = [0,y,0];
    y = 100*(y/nansum(y(:))); %normalize to % 
    patch(1:numel(y),y,fp.c_area(id(i),:),'edgecolor','none','facealpha',0.5)
    ylim([0 0.5]);
    set(gca,'clipping','off','ytick',[0 0.5])
    fp.FormatAxes(gca); box off
    xlabel('neurons');
    ylabel('% neurons');
    title(sprintf('%s',area_label{i}),'fontweight','normal');
    plot(get(gca,'xlim'),[0.25,0.25],'LineStyle',':','color',[0.6 0.6 0.6],'linewidth',1.2)
    xlim([1,numel(b)]) 
    camroll(-90)
    fp.FigureSizing(gcf,[3 2 0.75 4],[2 10 15 15])
    set(gca,'YAxisLocation','right')
end

end