function Figure9a(proj1,proj2)
%proj1 and 2 are the projections of activity during motif A and B,
%respectively onto the local region and subspace between regions

fp = fig_params_cortdynamics;
figure; hold on;
%plot the line of best fit
[~,xint] = AddLSline(proj1(:,1),proj1(:,2),proj1(:,1),[0.8 0.1 0.1],1.5);
AddLSline(proj2(:,1),proj2(:,2),proj2(:,1),[0.1 0.1 0.8],1.5);
plot(proj1(:,1),proj1(:,2),'o','color',[0.8 0.1 0.1],'markersize',fp.markersizesmall,'markerfacecolor',[0.8 0.1 0.1],'linestyle','none','linewidth',1.5);
plot(proj2(:,1),proj2(:,2),'o','color',[0.1 0.1 0.8],'markersize',fp.markersizesmall,'markerfacecolor',[0.1 0.1 0.8],'linestyle','none','linewidth',1.5);
%add a line at 45degrees 
fx = @(x) 1*x+xint;
unity = fx(proj1(:,1));
plot(proj1(:,1),unity,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle','-');
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 6 6],[10 10 12 12])

%add lines showing the distribution along each axis
xvals = get(gca,'xlim');
plot([xvals(1)+0.015,xvals(1)+0.015],[min(proj1(:,2)),max(proj1(:,2))],'linewidth',2,'color',[0.8 0.1 0.1])
plot([xvals(1)+0.025,xvals(1)+0.025],[min(proj2(:,2)),max(proj2(:,2))],'linewidth',2,'color',[0.1 0.1 0.8])
yvals = get(gca,'ylim');
plot([min(proj1(:,1)),max(proj1(:,1))],[yvals(1)+0.03,yvals(1)+0.03],'linewidth',2,'color',[0.8 0.1 0.1])
plot([min(proj2(:,1)),max(proj2(:,1))],[yvals(1)+0.05,yvals(1)+0.05],'linewidth',2,'color',[0.1 0.1 0.8])

xlabel({'Population Projection within','SS region (PC 1)'})
ylabel({'Subspace Projection between','SS and WHS (Dim 1)'})
title({'Example relationship between local','alignment and subspace activity'},'fontweight','normal')


end