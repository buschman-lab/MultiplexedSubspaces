function stats = Plot_CompareMotifStatistic(x,y,varargin)
%plotting options; 
opts.paired = 1; 
opts.yvals = [0 0.05];
opts.xvals = [0 0.05];
opts.pos = [3 3 4 5]; 
opts.figpos = [652   272   502   518];

opts = ParseOptionalInputs(opts,varargin);

stats = [];
figure; hold on; 
plot(linspace(opts.xvals(1),opts.xvals(2),1000),linspace(opts.xvals(1),opts.xvals(2),1000),'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--')
scatter(nanmean(x,1),nanmean(y,1),20,'filled','k')
%Add circles around the significant ones, were significance color
for i = 1:size(y,2)
    if opts.paired
        [pval, h] = signrank(x(:,i),y(:,i));
    else
        [pval, h] = ranksum(x(:,i),y(:,i));
    end
    stats.pval_store(i)=pval;
    if h == 1                
        if pval<0.0001
            scatter(nanmean(x(:,i)),nanmean(y(:,i)),100,[1 0 0],'filled')
        elseif pval<0.05%/size(y,2)
            scatter(nanmean(x(:,i)),nanmean(y(:,i)),50,[1 0 0],'filled')     
        end  
    end
    text([nanmean(x(:,i)),nanmean(x(:,i))],...
    [nanmean(y(:,i)),nanmean(y(:,i))],...
    sprintf('%d',i),'Fontsize',16,'FontWeight','normal',...
    'Color','k','HorizontalAlignment','center')  
end
ylim(opts.yvals);
xlim(opts.xvals);
setFigureDefaults;
% set(gca,'Position',opts.pos);
% set(gcf,'Position',[652   272   502   518])
