function PlotExamplePredictors(data)
fp = fig_params_cortdynamics;

area_label = data{1}(5).area_label;
x = cat(1,data{1}(5).area_val{1:7});
x = normalizeToBaseline(x,[1:2],'mean');
x = x(:,3:end,:);
x = x-nanmean(x,3);
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
area_label = area_label(1:7);
area_val = data{1}(5).area_val(1:7);

col = fp.c_area(1:7,:);

figure; hold on;
imagesc(x,[0 0.02]); 
set(gcf, 'Renderer', 'painters')
colormap(gca,flipud(gray));
xlim([-5,size(x,2)+0.5])
ylim([0,size(x,1)+0.5])

for i = 1:numel(area_label)    
    if i>1
        temp = cat(1,area_val{1:i-1});
        idx = [size(temp,1)+1,size(temp,1)+size(area_val{i},1)];
    else
        idx = [1,size(area_val{i},1)];
    end
    plot(idx,[-5 -5],'color',col(i,:),'linewidth',3); 
    text(idx(1)+(idx(2)-idx(1))/2, -10, area_label{i},'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end
xlim([0,size(x,2)+0.5])
ylim([-15,size(x,1)+0.5])
set(gca,'ytick',[],'xtick',[]);
fp.FormatAxes(gca); box on
set(gca,'Clipping','off')

fp.FigureSizing(gcf,[3 2 2 4],[2 10 30 15])


%plot example activity (average?)
figure; hold on; 
rng('default');
x = cat(1,area_val{1:7}); 
x = normalizeToBaseline(x,[1:2],'mean');
x = x(:,3:end,:);
x = x-nanmean(x,3);
x = squeeze(x(randperm(size(x,1),10),:,:));
for i = 1:5
   y = squeeze(x(i,:,:));
   y = abs(y);
   shadedErrorBar(1:12,nanmean(y,2),sem(y,2),'lineprops',{'color',[0.1 0.1 0.8 0.75],'linewidth',2});
end

%% the predicted

area_label = data{1}(5).area_label;
x = cat(1,data{1}(5).area_val{8});
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
area_label = area_label(8);
area_val = data{1}(5).area_val(8);

col = fp.c_area(8,:);

figure; hold on;
imagesc(x,[0 0.04]); 
set(gcf, 'Renderer', 'painters')
colormap(gca,flipud(gray));
xlim([-5,size(x,2)+0.5])
ylim([0,size(x,1)+0.5])

for i = 1:numel(area_label)    
    if i>1
        temp = cat(1,area_val{1:i-1});
        idx = [size(temp,1)+1,size(temp,1)+size(area_val{i},1)];
    else
        idx = [1,size(area_val{i},1)];
    end
    plot(idx,[-5 -5],'color',col(i,:),'linewidth',3); 
    text(idx(1)+(idx(2)-idx(1))/2, -10, area_label{i},'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','Color',fp.c_area(i,:),'fontsize',fp.font_size,'FontName',fp.font_name)
    
end
xlim([0,size(x,2)+0.5])
ylim([-15,size(x,1)+0.5])
set(gca,'ytick',[],'xtick',[]);
fp.FormatAxes(gca); box on
set(gca,'Clipping','off')

fp.FigureSizing(gcf,[3 2 1.5 4],[2 10 30 15])

end