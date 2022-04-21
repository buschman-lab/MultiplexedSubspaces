function [auc_all,pval_all] = SubspaceROC(data_rrr)
%Camden MacDowell - timeless

%loop through each target area
verbose=0;
%normalized PEV
n = numel(data_rrr(1).cvl_rrr);
auc_all = cell(1,n);
pval_all = cell(1,n);
for cur_a = 1:n
data = arrayfun(@(n) (1-data_rrr(n).cvl_rrr{cur_a})/max((1-data_rrr(n).cvl_rrr{cur_a}),[],'all'), 1:size(data_rrr,2),'UniformOutput',0);

%get the AUC of for each motif comparison and cross validation
auc = NaN(size(data,2),size(data,2),10);
pval = NaN(size(data,2),size(data,2));
for i = 1:size(data,2)
    for j = 1:size(data,2)
        if i~=j
           auc(i,j,:) = arrayfun(@(n) trapz([0,data{i}(n,:),1],[0,data{j}(n,:),1]),1:10);
           if nanmean(auc(i,j,:),'all')>=0.5
               pval(i,j) = signrank(squeeze(auc(i,j,:)),0.5,'tail','right');
           else
               pval(i,j) = signrank(squeeze(auc(i,j,:)),0.5,'tail','left');
           end
        else
           auc(i,j,:) = NaN(1,10);
           pval(i,j) = NaN;
        end
    end    
end

auc_all{cur_a} = auc;
pval_all{cur_a} = pval;

%visualize the mean
figure; hold on;
a = imagesc(nanmean(auc,3)); 
c=colorbar; colormap magma
ylabel(c,'AUC | light: y motif lower D than x motif')
p = pval;
p(p>0.05)=0.25;
p(p<0.25)=1;
p(isnan(p))=0;
a.AlphaData = p;
set(gca,'xlim',[0.5 size(data_rrr,2)+0.5],'ylim',[0.5 size(data_rrr,2)+0.5]);
title(['Area ',data_rrr(1).area_label{cur_a}],'FontWeight','normal');
xlabel('"x" motif'); ylabel('"y" motif');

if verbose
    %plot example curve
    t = tiledlayout(5,5);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    count = 1;
    for i = 1:size(data,2)
        if count>25; break; end                
        for j = 1:size(data,2)
            if count>25; break; end                
            if j>i
                nexttile
                hold on;
                arrayfun(@(n) plot([0,data{i}(n,:),1],[0,data{j}(n,:),1],'linewidth',0.25,'color',[0.75 0.75 0.75]),1:10)
                plot([0,nanmean(data{i}),1],[0,nanmean(data{j}),1],'linestyle','-','marker','none','linewidth',2,'color','k'); hold on;        
                xlabel(sprintf('Motif %d',i));
                ylabel(sprintf('Motif %d',j));
                plot(0:1,0:1,'linewidth',1,'color','r'); 
                if pval(i,j)<0.05
                    title(sprintf('auc=%0.2f',nanmean(auc(i,j,:))),'fontweight','bold')
                else
                    title(sprintf('auc=%0.2f',nanmean(auc(i,j,:))),'fontweight','normal')
                end
                axis square
                count = count+1;
            end        
        end
    end
    set(gcf,'units','normalized','position',[0 0 1 1])
end %verbose

end %brain region loop

end %function








% 
% figure; hold on; 
% plot(0:1,0:1)
% for i = 1:10
%    trapz([0,data{1}(i,:),1],[0,data{2}(i,:),1])
%    plot(data{1}(i,:),data{2}(i,:),'linewidth',1,'color','k','marker','o')
%    pause();
% end
% title('ROC')
% xlabel('Relative PEV Motif 1');
% ylabel('Relative PEV Motif 2');
% plot([0,nanmean(data{1}),1],[0,nanmean(data{2}),1],'linestyle','-','marker','o','markersize',3); hold on;
% trapz([0,nanmean(data{1}),1],[0,nanmean(data{2}),1]);
% %compare each motif to all other motifs
% a = data{1}; 
% b = data{2};
% temp = a-b;
% cvLoss = [mean(temp); std(temp)/sqrt(size(temp,1)) ];
% 
% p=[];
% for i = 1:30
%    p(i) = ranksum(a(:,i),b(:,i),'tail','left');
% end
% 
% figure; hold on;
% title('Difference in Relative PEV')
% ylabel('Motif 1 minus Motif 2 PEV');
% xlabel('Dimensions');
% y = cvLoss(1,:);
% e = cvLoss(2,:);
% errorbar(1:numel(y), y, e, 'o--','color',[0.8 0.1 0.1],'MarkerFaceColor','none','MarkerSize',8,'linewidth',1)
% plot([0 30],[0 0])
% 
% figure; hold on;
% title('ROC')
% xlabel('Relative PEV Motif 1');
% ylabel('Relative PEV Motif 2');
% plot(nanmean(a),nanmean(b),'linestyle','-','marker','o','markersize',3); hold on;
% plot(0:1,0:1)
%    
% temp = arrayfun(@(n) (1-data_rrr(n).cvl_rrr{1}), 1:size(data_rrr,2),'UniformOutput',0);
% a = temp{1}; 
% b = temp{2};
% figure; hold on;
% title('Absolute ROC')
% xlabel('Absolute PEV Motif 1');
% ylabel('Absolute PEV Motif 2');
% plot(nanmean(a),nanmean(b),'linestyle','-','marker','o','markersize',3); hold on;
% plot([0,0.25],[0,0.25])







