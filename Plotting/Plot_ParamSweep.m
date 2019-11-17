function fh = Plot_ParamSweep(in_fns,fh)

if nargin <2; fh = figure; end
set(0, 'currentfigure', fh);
hold on;

%load figure class
fp = fig_params; 

%load data
[ev_norm, param] = AnalyzeParamSweep(in_fns,{'L'});

%plot
shadedErrorBar(param*75,mean(ev_norm),sem(ev_norm,1),'lineProps',{'color',fp.c_discovery,...
    'linewidth',fp.dl_line_width},'patchSaturation',fp.dl_alpha);

%add line showing chosen length
line([13*75 13*75],[min(ev_norm(:)) 1],'linewidth',fp.dl_line_width,'linestyle','--',...
    'color',[0.5 0.5 0.5])

%add titles and format
xlabel('Motif Length (ms)');
ylabel({'Percent Explained Variance','(normalized)'});
fp.SetTitle(gca,{'Explained Varaince Plateues  '});
fp.FormatAxes(gca)

end