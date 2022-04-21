function RidgeRegressionCombinedFigures(data)
fp = fig_params_cortdynamics;

%motif_num are the two motifs to compare
[full_mdl,~] = LoadVariable(data,'ridge_performance',[]);

%plot performance across all areas and recordings relative to the
%trial-permuted distribution
figure; hold on; 
histogram(full_mdl(:),'BinWidth',0.025,'FaceAlpha',0.5,'FaceColor',fp.c_none,'EdgeAlpha',1,'EdgeColor',fp.c_none)
xlim([0 max(full_mdl(:))+0.025])
xlabel('Performance');
ylabel('Models');
fp.FormatAxes(gca); 
box on;
title({'trial-to-trial varaibility is','shared across brain areas'},'FontWeight','normal')
fp.FigureSizing(gcf,[5 2 2 2.5],[2 10 12 12])

end %function