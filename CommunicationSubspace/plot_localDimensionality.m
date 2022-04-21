function plot_localDimensionality(area_label,d,n,dimflag)
%camden - timeless
%see localDimensionality.m
if nargin <4; dimflag=1; end
figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on; 
subplot(2,1,1); 
bar(d); 
set(gca,'XTickLabel',area_label)
legend('location','northeastoutside')
if dimflag
    ylabel('optimal # of dimensions')
else
    ylabel('# of dimensions')
end
    
subplot(2,1,2); 
bar(n); 
set(gca,'XTickLabel',area_label)
legend('location','northeastoutside')
ylabel('# of neurons')

sgtitle('local dimensionality');

end %function end