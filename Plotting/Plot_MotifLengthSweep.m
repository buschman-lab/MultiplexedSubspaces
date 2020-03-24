function fh = Plot_MotifLengthSweep(in_fns,fh)

if nargin <2; fh = figure; end
set(0, 'currentfigure', fh);
hold on;

%load figure class
fp = fig_params; 

%load data
[ev_norm, param] = AnalyzeParamSweep(in_fns,{'L'});
ev_norm = ev_norm*100;

%plot
shadedErrorBar(param*75,mean(ev_norm),sem(ev_norm,1),'lineProps',{'color',[0.5 0.5 0.5],...
    'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);

%add line showing chosen length
line([13*75 13*75],[0 100],'linewidth',fp.dl_line_width,'linestyle','--',...
    'color',[0.1 0.1 0.1])

%add titles and format
xlabel('Motif Length (ms)');
ylabel({'Percent Explained','Variance (normalized)'});
xlim([75 5000]);
ylim([85 100]);
setFigureDefaults;
set(gca,'position',[3,5,6,6])
title({'~1s Motif Duration Reflects Plateau in';'Variance Captured by Motifs'},...
    'FontName','Arial','FontSize',16,'FontWeight','normal');
set(fh,'position',[680   400  450   600]);

end