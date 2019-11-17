
function pt = Plot_Snippet(data,x,col)
    y = convn(data,ones(6,1)/6,'same');
    y_mean = nanmean(y,2);
    y_sem = sem(y,2);
    shadedErrorBar(x,y_mean,y_sem,'lineprops',{'color',col});
    pt = plot(x,y_mean,'color',col);
end