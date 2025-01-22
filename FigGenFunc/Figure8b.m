function Figure8b(theta,ci)
fp = fig_params_cortdynamics;
figure; hold on;
histogram(theta,'facecolor',[0.5 0.5 0.5],'edgecolor','none','BinWidth',2.5)
yvals = get(gca,'ylim');
plot([0 0],yvals,'color','k','linewidth',1)
mu = rad2deg(circ_mean(deg2rad(theta)));
plot([mu mu],yvals,'color',[0.8 0.1 0.1],'linewidth',1)
plot([ci(1) ci(1)],yvals,'color',[0.8 0.1 0.1],'linewidth',1,'linestyle','--')
plot([ci(2) ci(2)],yvals,'color',[0.8 0.1 0.1],'linewidth',1,'linestyle','--')
xlabel('\Delta theta (degrees)');
ylabel('count');
title({'Magnitude of alignment','across recordings','(p<0.001)'},'FontWeight','normal')
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 3 4],[10 10 14 10])
set(gca,'ylim',yvals)

end