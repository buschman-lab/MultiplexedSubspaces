function Plot_CorticalNetworkSize(data)
%Camden MacDowell - timeless
fp = fig_params_cortdynamics;

%%
netsz = NaN(6,14,8,10);
for cur_rec = 1:6
    for cur_motif = 1:14
        for cur_a = 1:8         
            netsz(cur_rec,cur_motif,cur_a,:) = networkSize(data,cur_rec,cur_a,cur_motif);
        end
    end
end

%% plot and example for vis, motif 5 across recordings
figure; hold on; 
col = fp.c_area; col = col(8,:);
x = squeeze(netsz(:,5,8,:));
arrayfun(@(n) plot(1:size(x,2), x(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.5,'markersize',4), 1:size(x,1));
plot(1:size(x,2), nanmean(x),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');
xlim([0 size(x,2)]); 
xlabel('Subspace dimension');
ylabel({'% of Dorsal Cortex'});
title(sprintf('rec 1\n %s | motif %d','VIS',5),'fontweight','normal')
fp.FormatAxes(gca);  box on;
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])

%% Now do across all recordings and all motifs (per region)
close all;
area_label = data{1}(1).area_label;
x = squeeze(netsz);
ndim=10; 
col = fp.c_area;
figure; hold on;
for i = 1:size(col,1)   
    y = squeeze(x(:,:,i,:));
    y = permute(y,[3,2,1]);
    y = reshape(y,size(y,1),size(y,2)*size(y,3))';
    shadedErrorBar(1:ndim,nanmean(y),sem(y,1),'lineprops',{'color',[col(i,:),0.50],'linewidth',2});
    xlim([1 ndim]);
end
legend(area_label,'location','bestoutside');
xlabel('Subspace dimension');
ylabel({'% of Dorsal Cortex'});
title(sprintf('Cortical network Size | %s ',area_label{i}),'Fontweight','normal')
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 2 4 4],[10 10 15 10])

%% Plot the overlap across recordings for an example recording
[~,x,xx] = arrayfun(@(n) networkSize(data,n,8,5),1:6,'UniformOutput',0);

%replace nans with 0 where there are pixels
yy = cat(4,xx{:});
yy = sum(yy,4);

%for each dimension get the overlap
y = cat(4,x{:});
y(~isnan(y))=1;
y(isnan(y))=0;
y = sum(y,4);
y = double(y);
y(isnan(yy)) = NaN;

for i = 1:10
    figure; hold on;    
    PlotMesoFrame(y(:,:,i));
    set(gca,'CLim',[0 6])
    colorbar;
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    set(gca,'ydir','reverse');
    title(sprintf('motif %d | dim %d | area %s',motif, i, area));  
end

%% across dimensions get spatial correlation
cur_d=1:5;
rho = NaN(14,8,6,5);
rho_avg = NaN(14,8,6,5);
for cur_motif = 1:14
    for cur_area = 1:8
        [~,~,xx] = arrayfun(@(n) networkSize(data,n,cur_area,cur_motif),1:6,'UniformOutput',0);
        if cur_area == 4
            xx(4)=[];
            xx(3)=[];
        elseif cur_area == 6
            xx(3)=[];
        end
        
        xx = cellfun(@(x) reshape(x,68*68,10),xx,'UniformOutput',0);
        %for each dimension get the similarity
        for cur_rec = 1:numel(xx)
           y = xx{cur_rec};
           r = cellfun(@(x) triu(abs(corr(y(:,cur_d),x(:,cur_d),'rows','complete')),1),xx,'UniformOutput',0);
           r(cur_rec)=[];
           r = cat(3,r{:});
           r(r==0)=NaN;
           r = reshape(r,size(r,1)*size(r,2),size(r,3));
           rho(cur_motif,cur_area,cur_rec,1:size(r,2)) = nanmean(r);  
           
           %same thing across mismatch dimensions
           temp = cellfun(@(x) fliplr(x),xx,'UniformOutput',0);
           r = cellfun(@(x) triu(abs(corr(y(:,cur_d),x(:,cur_d),'rows','complete')),1),temp,'UniformOutput',0);
           r(cur_rec)=[];
           r = cat(3,r{:});
           r(r==0)=NaN;
           r = reshape(r,size(r,1)*size(r,2),size(r,3));
           rho_avg(cur_motif,cur_area,cur_rec,1:size(r,2)) = nanmean(r);
        end                    
    end
end

%% plot the correlation
r = permute(rho,[2,1,3,4]);
r = reshape(r,size(r,1),size(r,2),size(r,3)*size(r,4));
r = reshape(r,size(r,1),size(r,2)*size(r,3));

%reorganize by decreasing median
[~,idxsort] = sort(nanmedian(r,2),'descend');
col = fp.c_area; 
figure; hold on; 
boxplot(r(idxsort,:)','Notch','on')
a = get(get(gca,'children'),'children');
t = get(a,'tag'); 
idx = find(strcmp(t,'Box')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', col(i,:),'LineWidth',1.5); 
end
idx = find(strcmp(t,'Median')==1);
for i = 1:numel(idx)
    set(a(idx(i)), 'Color', 'k','LineWidth',1.5); 
end
xlabel('Area');
set(gca,'xticklabel',area_label(idxsort),'XTickLabelRotation',45)
ylabel('Spatial Similarity Across Recordings');
fp.FormatAxes(gca);  box on; grid on
fp.FigureSizing(gcf,[3 3 6 6],[10 10 10 10])

[p,~,stats] = kruskalwallis(r(idxsort,:)',[],'off');
title(sprintf('p=%0.3d',p),'fontweight','normal')
statresults = multcompare(stats,'Display','off');



end %function 

function [netsz,rho_sig,rho_all] = networkSize(data,cur_rec,area,motif)
netsz = NaN(1,10);
rho_sig=NaN(68,68,10);
rho_all=NaN(68,68,10);
try
    x = data{cur_rec};
    idx = find(cat(1,x.cur_motif)==motif & cat(1,x.cur_a)==area);
    rho_all = x(idx).rho_all; 
    sig_thresh = x(idx).sig_thresh; 
    
    rho_sig = rho_all;
    for i = 1:10
       temp = rho_all(:,:,i);
       temp(abs(temp)<sig_thresh(i))=NaN;
       rho_sig(:,:,i) = temp;
       x = temp;
       y = rho_all(:,:,i);
       netsz(i) = 100*sum(~isnan(x(:)))/sum(~isnan(y(:)));   
    end
catch
end

end %function eend
