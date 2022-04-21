function [fh,rho,pval,x,y] = Plot_CorrelateValuesBetweenRecordings(x,y,type,fp,tail,varargin)
%x and y are rec x motif matrices of values
%type is the way you want to compare things (i.e. normalization methods,
%controls, etc

opts.color_flag=1; %color_flag: 0 = all same color 1 = color by rec
opts.focus_motifs = []; %fill in specific motifs with special symbol to highlight
opts.xlabel = 'x';
opts.ylabel = 'y';
opts.corrtype = 'pearson';
opts.col = [0.75 0.75 0.75];
opts.rec = 1; %recording to plot
opts.addjitter=0; %for dimensionality
opts = ParseOptionalInputs(opts,varargin);

[n,m] = size(x);

switch type
    case 'combo' %plot all record
    case 'zscore'
        x = zscore(x,0,2);
        y = zscore(y,0,2);      
    case 'maxdiv'
        x = x./max(x,[],2);
        y = y./max(y,[],2);   
    case 'average'
        x = nanmean(x,1); 
        y = nanmean(y,1);
    case 'single_rec'
        x = x(opts.rec,:);
        y = y(opts.rec,:);
    otherwise
        error('unknown type of plot');
end

if opts.addjitter == 1 %x axis
    x = x+(rand(size(x,2),1)/2)'-0.5;
elseif opts.addjitter == 2 %y axis
    y = y+(rand(size(y,2),1)/2)'-0.5;
elseif opts.addjitter == 3 %both axes
    x = x+(rand(size(x,2),1)/2)'-0.5;
    y = y+(rand(size(y,2),1)/2)'-0.5;
end

figure; hold on;
if opts.color_flag==2
%     col = getColorPalet(n);
    col = arrayfun(@(n) linspace(n,0.1,size(x,1)),opts.col,'UniformOutput',0);
    col = cat(1,col{:})';
    arrayfun(@(n) plot(x(n,:),y(n,:),'marker','o','linestyle','none','markersize',fp.markersizesmall,'color',col(n,:)),1:size(x,1));
elseif opts.color_flag==1
    arrayfun(@(n) plot(x(n,:),y(n,:),'marker','o','linestyle','none','markersize',fp.markersizesmall,'color',opts.col(n,:)),1:size(x,1));
else
    col = opts.col;  
    plot(x(:),y(:),'marker','o','linestyle','none','markersize',fp.markersizesmall,'color',col)
end

[rho,pval] = corr(x(:),y(:),'type',opts.corrtype,'rows','complete','tail',tail);
h = plot(fitlm(x(:),y(:))); 
fp.FormatLSline(h);

if ~isempty(opts.focus_motifs)
    p = arrayfun(@(n) plot(x(:,opts.focus_motifs(n)),y(:,opts.focus_motifs(n)),'marker',fp.markers{n},...
        'linestyle','none','markersize',fp.markersizesmall+1,'markerfacecolor',col(1,:),'markeredgecolor',col(1,:)), 1:numel(opts.focus_motifs),'UniformOutput',0);
    temp = arrayfun(@(x) num2str(x),opts.focus_motifs,'UniformOutput',0);
    legend([p{:}],temp);
else
    hh = line(nan, nan, 'Color', 'none');
    legend(hh, {sprintf('r=%0.3f \n p=%0.3f',rho,pval)},'FontSize',fp.font_size,'FontName',fp.font_name,'Box','off','Location', 'best') 
end

xlabel(opts.xlabel,'Interpreter','tex');
ylabel(opts.ylabel,'Interpreter','tex');
fp.FormatAxes(gca); 
box on;
fp.FigureSizing(gcf,[3 2 5 6],[10 10 10 10])
fh = get(gcf);

end %function 
