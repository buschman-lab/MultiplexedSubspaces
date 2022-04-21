function PlotRegressionWeightsByArea(b,m,mout,area_lbl,area_name,fp,tails)
if nargin <7; tails = zeros(1,size(b,4)-1); end

x = squeeze(b(:,m,strcmp(area_lbl,area_name),strcmp(area_lbl,area_name)==0));
y = squeeze(b(:,mout,strcmp(area_lbl,area_name),strcmp(area_lbl,area_name)==0));
x(x==0)=NaN;
y(y==0)=NaN;

%reorder by decreasing values
[~,idx] = sort(nanmean(x),'descend');
x = x(:,idx); y = y(:,idx);
xlabels = area_lbl(strcmp(area_lbl,area_name)==0);
xlabels = xlabels(idx);

for i = 1:numel(tails)
    bar(i-0.5,nanmean(x(:,i)),'facecolor',[0.8 0.1 0.1],'facealpha',0.5,'EdgeAlpha',0,'barwidth',0.3)
    bar(i,nanmean(y(:,i)),'facecolor',[0.1 0.1 0.8],'facealpha',0.5,'EdgeAlpha',0,'barwidth',0.3)
    plot([i-0.5,i],[x(:,i),y(:,i)],'linewidth',1.5,'color',[0.75 0.75 0.75],'marker','o','markersize',5)
    if tails(i) == 0
        [~,p] = ttest(x(:,i),y(:,i),'tail','both');
    elseif tails(i) == 1
        [~,p] = ttest(x(:,i),y(:,i),'tail','left');
    else
        [~,p] = ttest(x(:,i),y(:,i),'tail','right');
    end
    text(i-0.25,max(cat(1,x(:,i),y(:,i))),sprintf('%0.3f',p));
end
set(gca,'xtick',(1:size(x,2))-0.25,'XTickLabel',xlabels,'XTickLabelRotation',45)

fp.FormatAxes(gca); 
box on;
fp.FigureSizing(gcf,[3 2 8 6],[10 10 12 10])


end %function 