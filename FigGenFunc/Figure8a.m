function Figure8a(theta,thetacon)

%% WHS --> RSP
figure;
ViolinAnglePlot(theta{1}(2,:),[0.8 0.1 0.1],1);
ViolinAnglePlot(thetacon{1}(2,:),[0.1 0.1 0.8],2);
set(gca,'xlim',[0.5 2.5])
set(gca,'xtick',[1,2]);

temp = theta{1}(2,:)-thetacon{1}(2,:);
if sum(temp)>0; p = sum(temp<=0)/numel(temp); else; p = sum(temp>=0)/numel(temp); end %we want one-sided for all of these, just written this way so can easily use for either)
title(sprintf('WHS-->RSP \np=%0.4f',p),'fontweight','normal')
set(gca,'xticklabel',{'A','B'})
xlabel('Motif')
ylim([60,90])



for i = 1:3
    figure;
    ViolinAnglePlot(theta{i}(1,:),[0.8 0.1 0.1],1);
    ViolinAnglePlot(thetacon{i}(1,:),[0.1 0.1 0.8],2);
    set(gca,'xlim',[0.5 2.5])
    temp = theta{i}(1,:)-thetacon{i}(1,:);
    if sum(temp)>0; p = sum(temp<=0)/numel(temp); else; p = sum(temp>=0)/numel(temp); end %we want one-sided for all of these, just written this way so can easily use for either)
    set(gca,'xticklabel',{'A','B'})
    xlabel('Motif')
   
    if i == 1 
        title(sprintf('WHS-->SS \np=%0.4f',p),'fontweight','normal')
    elseif i == 2
        ylim([40,90])
        title(sprintf('SS-->WHS \np=%0.4f',p),'fontweight','normal')
    else
        title(sprintf('RSP-->WHS \np=%0.4f',p),'fontweight','normal')
    end
end

end

function [mu,err] = ViolinAnglePlot(theta,col,xpos)
fp = fig_params_cortdynamics;
CompareViolins(theta,fp,'plotspread',0,'divfactor',1,'xpos',xpos,'plotaverage',0,'col',col,'distWidth',0.5);
mu = rad2deg(circ_mean(deg2rad(theta)));
line([xpos-0.075 xpos+0.075],[mu mu],'Color',col,'LineWidth',3)
fp.FormatAxes(gca)
fp.FigureSizing(gcf,[3 3 2.5 4],[10 10 14 10])
ylabel('Angle (degrees)');
err = rad2deg(circ_std(deg2rad(theta),[],[],2));
hold on; 
yval = get(gca,'ylim');
if yval(2)>90
    set(gca,'ylim',[yval(1),90]);
end
end