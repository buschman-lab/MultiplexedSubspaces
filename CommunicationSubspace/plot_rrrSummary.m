function plot_rrrSummary(score,cvl_rrr,cvl_fa,invalidFlag,strength)
%camden - timeless

figure; hold on; 

%rrr
cvLoss = [ mean(cvl_rrr); std(cvl_rrr)/sqrt(size(cvl_rrr,1)) ];
d = ModelSelect(cvLoss, 1:size(cvl_rrr,2));
y = 1-cvLoss(1,:);
e = cvLoss(2,:);
errorbar(1:numel(y), y, e, 'o--','color',[0.8 0.1 0.1],'MarkerFaceColor','none','MarkerSize',8,'linewidth',1)
errorbar(1:d, y(1:d), e(1:d), 'o--','color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',10,'linewidth',1.5)

%factor
cvLoss = [ mean(cvl_fa); std(cvl_fa)/sqrt(size(cvl_fa,1)) ];
% d = ModelSelect(cvLoss, 1:size(cvl_fa,2));
y = 1-cvLoss(1,:);
e = cvLoss(2,:);
errorbar(1:numel(y), y, e, 'o--','color',[0.1 0.1 0.8],'MarkerFaceColor','none','MarkerSize',8,'linewidth',1)
% errorbar(1:d, y(1:d), e(1:d), 'o--','color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',10,'linewidth',1.5)

%full 
cvLoss = [ mean(score); std(score)/sqrt(size(score,1)) ];
errorbar(0, cvLoss(1),cvLoss(2), 'marker','^','color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10,'linewidth',1.5)

if invalidFlag %if demarkated as a subspace
    title('no subspace','fontweight','normal')
else
    if strength ==1
        title('Strong Valid subspace','fontweight','normal')
    elseif strength ==2
        title('Medium Valid subspace','fontweight','normal')
    elseif strength ==3
        title('weak Valid subspace','fontweight','normal')
    end
end

xlabel('Number of predictive dimensions')
ylabel('Performance')

legend('rrr','rrr optimal','fa','full','Location', 'northeastoutside')

set(gca,'xlim',[-0.05 size(cvl_rrr,2)+0.25],'xtick',[1:size(cvl_rrr,2)])

end
