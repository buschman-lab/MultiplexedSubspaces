function [y,e] = plot_rrrFit(score,cvl_rrr)
%camden - timeless
%plot the rrrFit normalized to the explainable variance

figure; hold on; 

%rrr
cvLoss = [1-mean(cvl_rrr); std(cvl_rrr)/sqrt(size(cvl_rrr,1)) ];

%normalize to the explainable variance
cvLoss = cvLoss/nanmean(score); 
y = cvLoss(1,:);
d = find(y>0.9);
e = cvLoss(2,:);
errorbar(1:numel(y), y, e, 'o--','color',[0.8 0.1 0.1],'MarkerFaceColor','none','MarkerSize',8,'linewidth',1)
errorbar(1:d, y(1:d), e(1:d), 'o--','color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',10,'linewidth',1.5)

xlabel('Number of predictive dimensions')
ylabel('% of Explainable Variance')

legend('rrr','rrr optimal','Location', 'northeastoutside')

set(gca,'xlim',[-0.05 size(cvl_rrr,2)+0.25],'xtick',[1:4:size(cvl_rrr,2)])

end
