function PlotCraniotomyVariance(activity_in,activity_contra,n_edges,xval)
%Camden MacDowell - timeless
%plots the histograms of pixels values in each craniotomy
%versus the contralateral pixels
%Run "CompareROIVariance.mat" first

%Plot comparing hemispheres
figure; hold on;
[s_counts,edges] = histcounts(activity_contra(:),n_edges);
s_counts = 100*s_counts/sum(s_counts); %noramlize to the number of pixels 
[c_counts,~] = histcounts(activity_in(:),'BinEdges',edges);
c_counts = 100*c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,5),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,5),'linewidth',2,'color',[0.5 0.5 0.5])
xlim(xval)
[~,p,~,stats] = vartest2(activity_contra(:),activity_in(:));
title(sprintf('f=%0.2g p=%0.2g',stats.fstat,p),'fontweight','normal');

end