function Figure6a(y)
fp = fig_params_cortdynamics;
PlotMesoFrame(y);
set(gca,'Clim',[0,.4]);
cc = colorbar;
set(cc,'ytick',get(cc,'ylim'));
fp.FigureSizing(gcf,[3 2 6 6],[10 10 12 11])
title('Network Distribution','fontweight','normal');

end






