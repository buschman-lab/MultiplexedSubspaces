function [stats] = Plot_CortNetworksSummaryStats(data)
fp = fig_params_cortdynamics;


%% plot the size of the significant network over time
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

%% plot by area
figure; hold on;
col = fp.c_area;
for i = 1:numel(unig)
    xx = x(grp==unig(i),:);
    xx(xx==0)=NaN; %to avoid this effect being due to occurances of no significant pixels
    shadedErrorBar(1:size(xx,2),nanmean(xx),sem(xx,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
end
%significant decrease across subsequent dimensions
p = ones(1,size(x,2));
for i = 2:size(x,2)
    [~,p(i)] = ttest(x(:,i-1),x(:,i),'tail','right');    
end
p = p<(0.05/size(x,2));
for i = 1:size(x,2)
    if p(i)
        plot(i,nanmean(x(:,i)+std(x(:,i),[],1)),'*','color','k','markersize',fp.markersizesmall);
    end
end

ylabel({'Network size','(% of dorsal cortex)'})
xlabel('Dimension')
set(gca,'xtick',[1:3:10],'ytick',[0:25:100],'xlim',[1,10])
fp.FormatAxes(gca); box on; grid on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])

%% amount of overlap within a subspace, across dimensions [first for VIS area for 1 recording]
map = cat(2,data{:});
%go through each one
x = cell(1,size(map,2));
rho= cell(1,size(map,2));
ov = cell(1,size(map,2));
for j = 1:size(map,2)
    y = map(j).rho_all;
    p = map(j).sig_thresh;    
    ysig = y;
    for i = 1:10
       temp = ysig(:,:,i);
       temp(abs(temp)<p(i))=0;
       temp(abs(temp)>0)=1;
       ysig(:,:,i) = temp;
    end
    
    %flip rho_all so strongest is positive
    for i = 1:size(y,3)
        if nansum(y(:,:,i),'all') < nansum(-1*y(:,:,i),'all')
            y(:,:,i)=-1*y(:,:,i);
        end
    end    
    % get the overlap between all combineations
    ysig = reshape(ysig,size(ysig,1)*size(ysig,2),size(ysig,3));
    ytemp = reshape(y,size(y,1)*size(y,2),size(y,3));
    x{j} = NaN(10,10);
    rho{j} = NaN(10,10);
    temp = NaN(size(ysig,1),100);
    %get the overlap
    COUNT=1;
    for i = 2:size(ysig,2)
        for ii = 2:size(ysig,2)
            if ii>i
                x{j}(i,ii) = (2*sum((ysig(:,i)+ysig(:,ii))==2))/(sum(ysig(:,ii)==1)+sum(ysig(:,i)==1));
                rho{j}(i,ii) = corr(ytemp(:,i),ytemp(:,ii),'rows','complete');
                temp(:,COUNT) = (ysig(:,i)+ysig(:,ii))==2;
                COUNT = COUNT+1;
            end
        end
    end  
    ov{j} = nansum(temp,2);
end

%% just get the statistics
%average orthogonality
r = cat(3,rho{:});
r = reshape(r,size(r,1)*size(r,2),size(r,3));
r(isnan(r(:,1)),:)=[];
r = fisherZ(r);
clear stats;
[~,tempstats] = pairedBootstrap(r(:),@nanmean);
stats.rhomean = fisherInverse(tempstats.mean);
stats.rhoci = fisherInverse(tempstats.ci);

figure; hold on; 
rtemp = nanmean(r);
histogram(fisherInverse(rtemp),'binwidth',0.05,'Edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
mu = fisherInverse(nanmean(rtemp));
yvals = get(gca,'ylim');
plot([mu,mu],yvals,'linewidth',1.5,'color','r','linestyle','-')
ci = bootci(1000,@nanmean,rtemp);
ci = fisherZ(ci);
plot([ci(1),ci(1)],yvals,'linewidth',1.5,'color','r','linestyle',':')
plot([ci(2),ci(2)],yvals,'linewidth',1.5,'color','r','linestyle',':')
ylabel('# of cortical networks')
xlabel('Spatial Correlation')
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])

%%
r = cat(3,x{:});
r = reshape(r,size(r,1)*size(r,2),size(r,3));
r(isnan(r(:,1)),:)=[];
[bootovr,tempstats] = pairedBootstrap(r(:),@nanmean);
stats.overlapmean = fisherInverse(tempstats.mean);
stats.overlapmeanci = fisherInverse(tempstats.ci);

% Figure overlap
figure; hold on; 
histogram(100*nanmean(r),'binwidth',5,'Edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
ylabel('# of cortical networks')
mu = nanmean(r(:))*100;
yvals = get(gca,'ylim');
plot([mu,mu],yvals,'linewidth',1.5,'color','r','linestyle','-')
ci = bootci(1000,@nanmean,r(:))*100;
% plot([ci(1),ci(1)],yvals,'linewidth',1.5,'color','r','linestyle',':')
% plot([ci(2),ci(2)],yvals,'linewidth',1.5,'color','r','linestyle',':')
xlabel('% spatial overlap')
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])

%%

figure; hold on; 
y = nanmean(cat(2,ov{:}),2);
y = reshape(y,68,68);
y(isnan(y))=0;
PlotMesoFrame(y); %get color scaling
set(gca,'ydir','reverse')
set(gca,'clim',[0 8]); colorbar
% imagesc(y,[0 10]); colorbar
% set(gca,'ydir','reverse'); colormap magma;
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 2 4 4],[2 10 10 10])
xlim([0 68]); ylim([0 68]); axis off


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

