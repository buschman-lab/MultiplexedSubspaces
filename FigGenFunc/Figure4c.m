function Figure4c(ovrlap)
%% do full histogram and overlay lines for each brain area. 
fp = fig_params_cortdynamics;

binsize = 0.1;
edges = 0:binsize:1;

N = NaN(numel(edges)-1,8);
for i = 1:8
    x = squeeze(ovrlap(:,i,:,:));
    N(:,i) = histcounts(x(:),edges);
    N(:,i) = N(:,i)/sum(N(:,i));
end

figure; hold on; 
histogram(ovrlap(:),edges,'Normalization','probability','FaceAlpha',0.25,'FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.25)
col = fp.c_area;
for i = 1:8
    y = N(:,i);
    plot(edges(1:end-1)+binsize/2,y,'linewidth',1,'color',col(i,:),'linestyle','-')
end
ylabel('Probability');
xlabel('% Overlap');
fp.FormatAxes(gca); box on; grid off
fp.FigureSizing(gcf,[3 2 3 3],[2 10 10 10])

ymax = get(gca,'ylim');
xavg = nanmean(ovrlap(:));
plot([xavg,xavg],[0,ymax(2)],'linewidth',1,'color','r','linestyle','--');

title({'Similarity of','subspace networks'},'fontweight','normal');

end







