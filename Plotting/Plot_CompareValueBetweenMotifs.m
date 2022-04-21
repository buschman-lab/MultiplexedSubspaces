function [xx,fh] = Plot_CompareValueBetweenMotifs(x,m,fp,tail,col)
if nargin <5; col = [0.5 0.5 0.5]; end
bar(1:numel(m),nanmean(x(:,m)),'facecolor',col,'facealpha',0.4,'EdgeAlpha',0)
arrayfun(@(n) plot(1:numel(m),x(n,m),'linewidth',1.5,'color',col,'marker','o','markersize',fp.markersizesmall),1:size(x,1))

xx = x(:,m);

if numel(m)<5 %don't draw crazy number of ttests
    a = nchoosek(m,2);
    b = nchoosek(1:numel(m),2);
    for i = 1:size(a,1)   
        [~,p] = ttest(x(:,a(i,1)),x(:,a(i,2)),'tail',tail);
        yval = max(x(:,m),[],'all');
        yval = yval+(i*0.05*yval);
        plot([b(i,1),b(i,2)],[yval,yval],'linewidth',1,'color','k')
        text(b(i,1)+abs(b(i,1)-b(i,2))/2,yval,sprintf('%0.2f',p),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
else %just do an anova
    [pval,tbl] = anova1(x,[],'off');
    hh = line(nan, nan, 'Color', 'none');
    legend(hh, {sprintf('fstat=%0.3f \n p=%0.23',tbl{2,5},pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best')     
end
xlabel('motif');
set(gca,'xtick',1:numel(m),'xticklabel',m); 
set(gca,'xlim',[0.25,numel(m)+0.75]);

fp.FormatAxes(gca); 
box on;
if numel(m)>4
    fp.FigureSizing(gcf,[3 2 numel(m)*1.5 4],[10 10 numel(m)*2 8])
else
    fp.FigureSizing(gcf,[3 2 numel(m)*1.5 4],[10 10 14 8])
end

fh = get(gcf);

end %function end