function Figure1f(mdl_perf)
fp = fig_params_cortdynamics;
%plot performance across all areas and recordings
figure; hold on; 
histogram(mdl_perf(:),'BinWidth',0.025,'FaceAlpha',0.5,'FaceColor',fp.c_none,'EdgeAlpha',1,'EdgeColor',fp.c_none)
xlim([0 max(mdl_perf(:))+0.025])
xlabel('Performance');
ylabel('Datasets');
fp.FormatAxes(gca); 
box on;
title({'Model Fit'},'FontWeight','normal')
fp.FigureSizing(gcf,[5 2 2 2.5],[2 10 12 12])

end %function