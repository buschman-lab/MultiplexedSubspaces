function [stats] = Plot_CortNetworksSummaryStats(data)
fp = fig_params_cortdynamics;




%% plot the size of the significant network over dimensions
map = cat(2,data{:});
%go through each one
x = NaN(size(map,2),10);
for j = 1:size(map,2)
    y = map(j).rho_all;
    p = map(j).sig_thresh;    
    p = map(j).sig_thresh;    
    for i = 1:10
        temp = abs(y(:,:,i));
        x(j,i) = 100*(sum(temp>p(i),'all')/sum(~isnan(temp(:))));        
    end
end

grp = cat(1,map(:).cur_a);
unig = unique(grp);

%% show that the first dimension spans the cortex
figure; hold on;
a = x(:,1);
b = x(:,2:10);
b(b<1)=[];
b = b(:);
bar(1,nanmean(a(:)),'FaceColor',fp.c_ff,'FaceAlpha',0.5,'edgecolor','none')
bar(2,nanmean(b(:)),'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.5,'edgecolor','none')
y = rand(1,size(a,1),1)/4+1-0.125;
plot(y,a(:),'marker','.','markersize',0.5,'linestyle','none','color',fp.c_ff)
y = rand(1,size(b,1),1)/4+2-0.125;
plot(y,b(:),'marker','.','markersize',0.5,'linestyle','none','color',[0.2 0.2 0.2])
% CompareViolins(a',fp,'plotspread',0,'divfactor',0.5,'plotaverage',1,'col',{fp.c_ff},'distWidth',0.75,'xpos',1);
% CompareViolins(b(:)',fp,'plotspread',0,'divfactor',0.5,'plotaverage',1,'distWidth',0.75,'xpos',2);
xlim([0.5 2.5])
ylim([0 100])
ylabel({'Network size','(% of dorsal cortex)'})
xlabel('Dimension')
set(gca,'xtick',[1,2],'ytick',[0:25:100])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])

% bootstrap to get their difference
rng('default');
d = NaN(1,1000);
for i = 1:1000    
    idx = datasample(1:numel(a),numel(a));
    idxb = datasample(1:numel(b),numel(b));
    d(i) = nanmean(a(idx))-nanmean(b(idxb));
end
p = sum(d<0)/numel(d);

title(sprintf('%0.3f',p))

%normalized histogram
figure; hold on;
a = x(:,1);
b = x(:,2:10);
b(b<1)=[];
b = b(:);
histogram(a,'binwidth',2,'Facecolor',fp.c_ff,'edgealpha',0,'Normalization','probability');
histogram(b,'binwidth',2,'Facecolor',[0.5 0.5 0.5],'edgealpha',0,'Normalization','probability');
xlabel({'Network size','(% of dorsal cortex)'})
ylabel('Fraction of Networks')
set(gca,'xtick',[0:25:100]);
legend('Dimension 1','Dimensions 2-10','Location','best');
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
yvals = get(gca,'ylim');
plot([nanmean(b(:)),nanmean(b(:))],yvals,'linestyle',':','color','k','linewidth',1.5)
plot([nanmean(a(:)),nanmean(a(:))],yvals,'linestyle',':','color',fp.c_ff,'linewidth',1.5)
ylim([yvals]);
[~,stats] = pairedBootstrap(b(:),@nanmean);
[~,stats2] = pairedBootstrap(a(:),@nanmean);
title(sprintf('%0.3f %0.3f %0.3f \n %0.3f %0.3f %0.3f',stats.mean,stats.ci(1),stats.ci(2),stats2.mean,stats2.ci(1),stats2.ci(2)),'fontweight','norma');

figure; hold on;
a = x(:,1);
b = x(:,2:10);
b(b<1)=[];
b = b(:);
a = cat(1,a,NaN(numel(b)-numel(a),1));
temp = [a,b]';
CompareViolins(temp,fp,'plotspread',0,'divfactor',1,'plotaverage',1,'col',{fp.c_ff,[0.2 0.2 0.2]},'distWidth',0.75,'sidebyside',1);
xlim([0.75 1.25])
ylim([0 100])
ylabel({'Network size','(% of dorsal cortex)'})
xlabel('Dimension')
set(gca,'xtick','','ytick',[0:25:100])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])

% bootstrap to get their difference
rng('default');
d = NaN(1,1000);
for i = 1:1000    
    idx = datasample(1:numel(a),numel(a));
    idxb = datasample(1:numel(b),numel(b));
    d(i) = nanmean(a(idx))-nanmean(b(idxb));
end
p = sum(d<0)/numel(d);

title(sprintf('%0.3f',p))

%%
% %% Compare the size of subcortical to cortical 
% map = cat(2,data{:});
% %go through each one
% x = cell(1,8);
% for cur_a = 1:8
%     idx = find(cat(1,map(:).cur_a)==cur_a);
%     tempmap = map(idx); 
%     for j = 1:size(tempmap,2)
%     y = tempmap(j).rho_all;
%     p = tempmap(j).sig_thresh;    
%     p = tempmap(j).sig_thresh;   
%     for i = 1:10
%         temp = abs(y(:,:,i));
%         x{cur_a}(j,i) = 100*(sum(temp>p(i),'all')/sum(~isnan(temp(:))));        
%     end
%     end
% end
% 
% figure; hold on
% for i = 1:8
%     xboot = x{i}(:);
%     xboot = pairedBootstrap(xboot,@nanmean);
%     CompareViolins(xboot',fp,'plotspread',0,'divfactor',0.2,'plotaverage',1,'distWidth',0.75,'xpos',i);    
% end
% xlim([0.5 8.5])
% 
% 
% 
% %% plot by area
% figure; hold on;
% col = fp.c_area;
% for i = 1:numel(unig)
%     xx = x(grp==unig(i),:);
%     xx(xx==0)=NaN; %to avoid this effect being due to occurances of no significant pixels
%     shadedErrorBar(1:size(xx,2),nanmean(xx),sem(xx,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
% end
% %significant decrease across subsequent dimensions
% p = ones(1,size(x,2));
% for i = 2:size(x,2)
%     [~,p(i)] = ttest(x(:,i-1),x(:,i),'tail','right');    
% end
% p = p<(0.05/size(x,2));
% for i = 1:size(x,2)
%     if p(i)
%         plot(i,nanmean(x(:,i)+std(x(:,i),[],1)),'*','color','k','markersize',fp.markersizesmall);
%     end
% end
% 
% ylabel({'Network size','(% of dorsal cortex)'})
% xlabel('Dimension')
% set(gca,'xtick',[1:3:10],'ytick',[0:25:100],'xlim',[1,10])
% fp.FormatAxes(gca); box on; grid on; 
% fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
% 
% %% amount of overlap within a subspace, across dimensions [first for VIS area for 1 recording]
% map = cat(2,data{:});
% %go through each one
% x = cell(1,size(map,2));
% rho= cell(1,size(map,2));
% ov = cell(1,size(map,2));
% for j = 1:size(map,2)
%     y = map(j).rho_all;
%     p = map(j).sig_thresh;    
%     ysig = y;
%     for i = 1:10
%        temp = ysig(:,:,i);
%        temp(abs(temp)<p(i))=0;
%        temp(abs(temp)>0)=1;
%        ysig(:,:,i) = temp;
%     end
%     
%     %flip rho_all so strongest is positive
%     for i = 1:size(y,3)
%         if nansum(y(:,:,i),'all') < nansum(-1*y(:,:,i),'all')
%             y(:,:,i)=-1*y(:,:,i);
%         end
%     end    
%     % get the overlap between all combineations
%     ysig = reshape(ysig,size(ysig,1)*size(ysig,2),size(ysig,3));
%     ytemp = reshape(y,size(y,1)*size(y,2),size(y,3));
%     x{j} = NaN(10,10);
%     rho{j} = NaN(10,10);
%     temp = NaN(size(ysig,1),100);
%     %get the overlap
%     COUNT=1;
%     for i = 1:size(ysig,2)
%         for ii = 1:size(ysig,2)
%             if ii>i
%                 x{j}(i,ii) = (2*sum((ysig(:,i)+ysig(:,ii))==2))/(sum(ysig(:,ii)==1)+sum(ysig(:,i)==1));
%                 rho{j}(i,ii) = corr(ytemp(:,i),ytemp(:,ii),'rows','complete');
%                 temp(:,COUNT) = (ysig(:,i)+ysig(:,ii))==2;
%                 COUNT = COUNT+1;
%             end
%         end
%     end  
%     ov{j} = nansum(temp,2);
% end
% 
% %% just get the statistics
% %average orthogonality
% r = cat(3,rho{:});
% r = reshape(r,size(r,1)*size(r,2),size(r,3));
% r(isnan(r(:,1)),:)=[];
% r = fisherZ(r);
% clear stats;
% [~,tempstats] = pairedBootstrap(r(:),@nanmean);
% stats.rhomean = fisherInverse(tempstats.mean);
% stats.rhoci = fisherInverse(tempstats.ci);
% 
% figure; hold on; 
% rtemp = nanmean(r);
% histogram(fisherInverse(rtemp),'binwidth',0.05,'Edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
% mu = fisherInverse(nanmean(rtemp));
% yvals = get(gca,'ylim');
% plot([mu,mu],yvals,'linewidth',1.5,'color','r','linestyle','-')
% ci = bootci(1000,@nanmean,rtemp);
% ci = fisherZ(ci);
% plot([ci(1),ci(1)],yvals,'linewidth',1.5,'color','r','linestyle',':')
% plot([ci(2),ci(2)],yvals,'linewidth',1.5,'color','r','linestyle',':')
% ylabel('# of cortical networks')
% xlabel('Spatial Correlation')
% fp.FormatAxes(gca); box on; 
% fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
% 
% %%
% r = cat(3,x{:});
% r = reshape(r,size(r,1)*size(r,2),size(r,3));
% r(isnan(r(:,1)),:)=[];
% [bootovr,tempstats] = pairedBootstrap(r(:),@nanmean);
% stats.overlapmean = fisherInverse(tempstats.mean);
% stats.overlapmeanci = fisherInverse(tempstats.ci);
% 
% % Figure overlap
% figure; hold on; 
% histogram(100*nanmean(r),'binwidth',5,'Edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
% ylabel('# of cortical networks')
% mu = nanmean(r(:))*100;
% sem(r(:))*100
% yvals = get(gca,'ylim');
% plot([mu,mu],yvals,'linewidth',1.5,'color','r','linestyle','-')
% ci = bootci(1000,@nanmean,r(:))*100;
% % plot([ci(1),ci(1)],yvals,'linewidth',1.5,'color','r','linestyle',':')
% % plot([ci(2),ci(2)],yvals,'linewidth',1.5,'color','r','linestyle',':')
% xlabel('% spatial overlap')
% fp.FormatAxes(gca); box on; 
% fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
% 
% %%
% 
% figure; hold on; 
% y = nanmean(cat(2,ov{:}),2);
% y = reshape(y,68,68);
% y(isnan(y))=0;
% PlotMesoFrame(y); %get color scaling
% set(gca,'ydir','reverse')
% set(gca,'clim',[0 13]); colorbar
% % imagesc(y,[0 10]); colorbar
% % set(gca,'ydir','reverse'); colormap magma;
% fp.FormatAxes(gca); box on; 
% fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
% xlim([0 68]); ylim([0 68]); axis off


%% Get the map of where things overlap
% close all;
% figure; hold on; 
% t=tiledlayout(2,4,'Padding','normal','TileSpacing','normal');
% %get mean across areas 
% for i =1:8
%     nexttile; hold on; 
%     idx = cat(1,map(:).cur_a)==i;    
%     y = reshape(nanmean(cat(2,ov{idx}),2),68,68);
%     imagesc(y); colorbar
%     set(gca,'ydir','reverse'); colormap magma;
%     title(sprintf('areas %d',i),'fontweight','normal');
% end





%%

