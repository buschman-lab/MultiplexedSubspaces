function fh = Plot_CompareBehaviorStates_2States(weight,yvals)

rng('default');
sig_color = [0.75 0.75 0.75; 0.25 0.25 0.25];
pval = [];
y=[];
y_diff = [];
for i = 1:size(weight,1)
   temp = cat(1,weight{i,:})';    
   %temp is a n x 2 matrix of average motif weights where columns are state 1 and 2 and
   %rows are instances
   pval(i) = ranksum(temp(:,1),temp(:,2));   
   temp = temp/nanmean(temp(:)); %get normalized differences  
   y(:,i) = nanmean(temp);
end

figure('Position',[0 0 1000 5000]); hold on; 
s1 = subplot(312,'Units','centimeters','Position',[8 6 10 2.0]); hold on

axes(s1);
imagesc(y,[yvals])
colormap(gca,flipud(redgreencmap(256,'Interpolation','quadratic')));
c = colorbar;
set(c,'YTick',[yvals(1), 1, yvals(2)]);
ylabel(c,{'Normalized Motif';'Intensity'},'FontSize',16,'Fontweight','normal','FontName','Arial');
set(c,'units','centimeters','position',[18.25 6 0.5 2])

for x_grid = 0.5:1:size(y,2)+0.5
    line([x_grid,x_grid],[0.5,size(y,1)+0.5],'linewidth',1.5,'color','w')    
end
for y_grid = 0.5:1:size(y,1)+0.5
    line([0.5,size(y,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
end  
xlabel('Basis Motifs')
%flip axes


xlim([0.5 size(y,2)+.5])
ylim([0.5 size(y,1)+0.5])
set(gca,'YColor','k')
setFigureDefaults

for i = 1:numel(pval)
   th = AddSig(1,pval(i),[i-0.1,i-0.1,2.6,2.6],1,0,1,90);
   set(th,'HorizontalAlignment','left')
%    if pval(i)<=0.05/numel(pval)
%        set(th,'fontweight','bold')
%    end
end

set(gcf,'Position',[680   150  875  650]);

fh = gcf;
end


%%
% 
% rng('default');
% sig_color = [0.75 0.75 0.75; 0.25 0.25 0.25];
% pval = [];
% y=[];
% y_diff = [];
% for i = 1:size(weight,1)
%    temp = cat(1,weight{i,:})'; 
%    
%    %temp is a n x 2 matrix of average motif weights where columns are state 1 and 2 and
%    %rows are instances
%    pval(i) = kruskalwallis(temp,[],'off');   
%    temp = temp/nanmean(temp(:)); %get normalized differences    
%    y(:,i) = nanmean(temp);
%    
%    %bootstrap distribution
%    boot_dist = NaN(1,1000);
%    for j = 1:5000 
%       temp_x = randsample(temp(:,1),numel(temp(:,1)),'true');
%       temp_y = randsample(temp(:,2),numel(temp(:,2)),'true');
%       boot_dist(j)=nanmean(temp_x)-nanmean(temp_y);
%    end
%    %get confidence interval and mean;
%    y_diff(1,i) = nanmean(boot_dist);
%    ci = [prctile(boot_dist,2.5),prctile(boot_dist,97.5)];
% 
%    
%    %SEM
%    y_diff(3,i) = y_diff(1,i)-ci(1);
%    y_diff(3,i) = ci(2)-y_diff(1,i); 
% 
% %    y_diff(1,i) = nanmean(temp(:,1)-temp(:,2));
% %    ci = bootci(1000,@nanmean,(temp(:,1)-temp(:,2)));
% %    y_diff(2,i) = y_diff(1,i)-ci(1);
% %    y_diff(3,i) = ci(2)-y_diff(1,i);   
% %    y_diff(2,i) = y_diff(1,i)-sem(temp(:,1)-temp(:,2));
% %    y_diff(3,i) = sem(temp(:,1)-temp(:,2))-y_diff(1,i);   
% end
% 
% figure('Position',[0 0 1000 5000]); hold on; 
% s1 = subplot(312,'Units','centimeters','Position',[8 6 10 2.0]); hold on
% s2 = subplot(311,'Units','centimeters','Position',[8 9 10 2.0]); hold on
% 
% axes(s1);
% imagesc(y,[yvals])
% colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));
% c = colorbar;
% set(c,'YTick',[yvals(1), 1, yvals(2)]);
% ylabel(c,{'Normalized Motif';'Intensity'},'FontSize',16,'Fontweight','normal','FontName','Arial');
% set(c,'units','centimeters','position',[18.25 6 0.5 2.5])
% 
% for x_grid = 0.5:1:size(y,2)+0.5
%     line([x_grid,x_grid],[0.5,size(y,1)+0.5],'linewidth',1.5,'color','w')    
% end
% for y_grid = 0.5:1:size(y,1)+0.5
%     line([0.5,size(y,2)+0.5],[y_grid, y_grid],'linewidth',1.5,'color','w')    
% end  
% xlabel('Basis Motifs')
% 
% 
% xlim([0.5 size(y,2)+.5])
% ylim([0.5 size(y,1)+0.5])
% set(gca,'YColor','k')
% setFigureDefaults
% 
% axes(s2); hold on
% cla
% %Plot the difference and the significance
% line([0.5 size(weight,1)+.5],[0 0],'linestyle','--','linewidth',1.5,'color',[0.1 0.1 0.1]);
% plot(y_diff(1,:),'.','MarkerSize',5,'MarkerEdgeColor',[0.4 0.4 0.4],'Color',[0.4 0.4 0.4])
% for i = 1:numel(pval)
%     errorbar(i,y_diff(1,i),y_diff(2,i),y_diff(3,i),'LineStyle','none','Linewidth',1.5,'Marker','.','MarkerSize',20,...
%         'Color',sig_color((pval(i)<=0.05/numel(pval))+1,:))
% end
% for i = 1:numel(pval)
%    AddSig(1,pval(i),[i-0.1,i-0.1,y_diff(1,i),y_diff(1,i)],1,4,1,90)
% end
% 
% %Change marker color for significant motifs (lighter)
% xlim([0.5 size(weight,1)+.5])
% set(gca,'XTick',(1:2:size(weight,2)),'Ytick',[min(get(gca,'Ytick')),max(get(gca,'Ytick'))]);
% ylabel({'\Delta\mu';''},'Rotation',0,'Units','Centimeters','position',[11.75 1.4]);
% set(gca,'yaxislocation','right')
% set(gca,'XColor','none');
% setFigureDefaults
% 
% 
% set(gcf,'Position',[680   150  875  650]);
% 
% fh = gcf;