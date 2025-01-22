function [m,xint, xval] = AddLSline(x,y,xval,c,w)

% p = polyfit(x,y,1); 
% f = polyval(p,xval); 
% plot(xval,f,'color',c,'linewidth',w);

xval = sort(xval,'ascend');
lm = fitlm(x,y);
[ypred,yci] = predict(lm,xval);
ci = abs(ypred-yci);
shadedErrorBar(xval,ypred,ci,'lineprops',{'color',[c 0.75],'linewidth',1});

%get the slope
lm=fitlm(x,y);             % a linear model
m=lm.Coefficients.Estimate;
xint = m(1);
m = m(2);

% p = plot(lm);
% delete(p(1));
% for i = 2:4
%     p(i).LineWidth=w;
%     p(i).Color=c;
% end



end