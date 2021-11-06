function PlotHemisphereVariance(dff,type,n_edges,xval)
%Camden MacDowell - timeless
%plots the histograms of pixels values in each hemisphere
%type can be pixelwise variance ('pxlvar') or all values over time ('all')

%mask is not perfectly centered (which is deliberate for earlier preprocessing legacy), but it means we have to shift it on this end.
R_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
R_hemi(:,37:end,:) = fliplr(dff(:,3:34,:));
R_hemi = reshape(R_hemi,[size(R_hemi,1)^2,size(R_hemi,3)]);
L_hemi = NaN(size(dff,1),size(dff,2),size(dff,3)); %convert back to full size
L_hemi(:,37:end,:) = dff(:,37:end,:);
L_hemi = reshape(L_hemi,[size(L_hemi,1)^2,size(L_hemi,3)]);

switch type
    case 'all'
        %do nothing
    case 'pxlstd'
        R_hemi = nanstd(R_hemi,[],2);
        L_hemi = nanstd(L_hemi ,[],2);
    case 'pxlvar'
        R_hemi = nanvar(R_hemi,[],2);
        L_hemi = nanvar(L_hemi ,[],2);               
    otherwise
        error('unknown type');
end

%Plot comparing hemispheres
figure; hold on;
[s_counts,edges] = histcounts(R_hemi(:),n_edges);
s_counts = s_counts/sum(s_counts); %noramlize to the number of pixels 
[c_counts,~] = histcounts(L_hemi(:),'BinEdges',edges);
c_counts = c_counts/sum(c_counts);
edges = edges(1:end-1)+diff(edges)/2;
bar(edges,c_counts,'FaceColor',[0.9100    0.4100    0.1700],'FaceAlpha',0.5,'EdgeAlpha',0); 
bar(edges,s_counts,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0); 
plot(edges,smooth(c_counts,10),'linewidth',2,'color',[0.9100    0.4100    0.1700])
plot(edges,smooth(s_counts,10),'linewidth',2,'color',[0.5 0.5 0.5])
xlim([xval])
[~,p,~,stats] = vartest2(R_hemi(:),L_hemi(:));
title(sprintf('f=%0.2g p=%0.2g',stats.fstat,p),'fontweight','normal');

end