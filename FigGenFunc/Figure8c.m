function Figure8c(theta,xpos)
fp = fig_params_cortdynamics;
figure; hold on;
for i = 1:numel(theta)    
    if i==1
        CompareViolins(theta{i}',fp,'plotspread',0,'divfactor',0.4,'xpos',xpos(i),'col',{fp.c_ff},'distWidth',1.4);    
    else
        CompareViolins(theta{i}',fp,'plotspread',0,'divfactor',0.55,'xpos',xpos(i),'col',{[0.25 0.25 0.25]},'distWidth',0.9);%0.5
    end   
end
camroll(-90)
set(gca,'xtick',xpos)
xlim([0 xpos(end)+0.75])
plot([0 xpos(end)+1],[0 0],'linestyle','-','color',[0.3 0.3 0.3],'linewidth',1.5)
ylim([-10 5]);
ylabel('\Delta theta');
set(gca,'YAxisLocation','right')
set(gca, 'YGrid', 'off', 'XGrid', 'on','GridAlpha',0.15);
fp.FormatAxes(gca); box on; 
fp.FigureSizing(gcf,[3 3 4 5.75],[10 10 14 10])
set(gca,'XTickLabel',{'All','HPC,VIS,WHS','VIS,PL,WHS','SS,WHS,RSP','PL,WHS,FMR','FMR,SS,VIS','TH,WHS,FMR'})
end
