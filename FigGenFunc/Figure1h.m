function Figure1h(rrr_mdl, area_label)
fp = fig_params_cortdynamics;

ndim = 15;
thresh = 0.80;

%plot the relative activity 
col = fp.c_area; col = col(strcmp(area_label,'VIS'),:);
figure; hold on;
plot([0 ndim],[thresh thresh] ,'linestyle','--','color','k','linewidth',1.5);
rrr_mdl = reshape(rrr_mdl(:,:,1:ndim),size(rrr_mdl,1)*size(rrr_mdl,2),ndim);
arrayfun(@(n) plot(1:size(rrr_mdl,2), rrr_mdl(n,:),'linestyle','-','marker','none','color',[0.75 0.75 0.75],'linewidth',0.5,'markersize',4), 1:size(rrr_mdl,1));

plot(1:size(rrr_mdl,2), nanmean(rrr_mdl),'linestyle','-','marker','o','color',col,'linewidth',1.5,'markersize',5,'MarkerFaceColor','none');  
xlim([0 ndim]); ylim([0 1]);
xlabel('# of dimensions');
ylabel({'Performance','(r^2 of explainable variance)'});
title(sprintf('%s ','Visual Region'),'fontweight','normal')
fp.FormatAxes(gca);  box on;
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])

end