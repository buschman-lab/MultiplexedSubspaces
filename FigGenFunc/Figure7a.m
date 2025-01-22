function Figure7a(a,b)
% areaname = {'WHS','SS','RSP'};
fp = fig_params_cortdynamics;
%get all motifs
rng('default');
for i = 1:3
    %simple mean + sem plot
    figure; hold on; 
    plot([1,2],[nanmean(a{i}),nanmean(b{i})],'color',[0.5 0.5 0.5],'linewidth',1)
    errorbar(1,nanmean(a{i}),sem(a{i}),'marker','o','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall);
    errorbar(2,nanmean(b{i}),sem(b{i}),'marker','o','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall);
    
    p = signrank(a{i},b{i}); 
    set(gca,'XTickLabel',{'A','B'},'xlim',[0.5 2.5])
    xlabel('Motif')
    yval = get(gca,'ylim');
    set(gca,'ylim',[round(yval(1)-0.025,2),round(yval(2)+0.025,2)]);
    fp.FormatAxes(gca);  box on; grid on; axis square; 
    fp.FigureSizing(gcf,[3 3 1.25 2],[10 10 14 10])   
    yval = get(gca,'ylim');
    set(gca,'ytick',yval);
    ylabel({'Neural','Activity'});    

end

end




