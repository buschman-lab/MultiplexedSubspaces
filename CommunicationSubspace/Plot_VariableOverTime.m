function y = Plot_VariableOverTime(y,c,type)
%y is rec x time matrix
if nargin <3; type = 'max'; end

if size(y,2)==1
    y = y';
end

x = 1:size(y,2);
switch type
    case 'max'        
        y = y./nanmax(y,[],2);
    case 'minmax'
        y = (y-nanmin(y,[],2))./(nanmax(y,[],2)-nanmin(y,[],2)); 
    case 'minmaxminusmean'
        y = (y-nanmin(y,[],2))./(nanmax(y,[],2)-nanmin(y,[],2)); 
        y = y-nanmean(y,2);
    case 'zscore'
        y = (y-nanmean(y,2))./nanstd(y,[],2);
    case 'none'
        
    otherwise 
        error('no type');
end

if size(y,1)>1
    shadedErrorBar(x,nanmean(y),sem(y),'lineprops',{'color',c,'linewidth',1.5});
else
    plot(x,y,'color',c,'linewidth',1.5);
end

xlabel('time (ms)'); set(gca,'Xtick',0:2:10,'XTickLabel',floor([0:2:10]*133));

end