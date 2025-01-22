function Figure7b(a,b,area)
   fp = fig_params_cortdynamics;
   
   i=1; %PC1
   %plot planes
   x = cat(3,a{:});
   y = cat(3,b{:});
   x = squeeze(x(:,i,:));   
   y = squeeze(y(:,i,:));
   figure; hold on; 
   x = movmean(x,2,1);
   y = movmean(y,2,1);
   [X,Y,Z] = Get3DPlane(x,0.75);
   patch(X,Y,Z,[0.8 0.1 0.1],'EdgeColor','none','FaceAlpha',.5);
   [X,Y,Z] = Get3DPlane(y,0.75);
   patch(X,Y,Z,[0.1 0.1 0.8],'EdgeColor','none','FaceAlpha',.5);   
   xlabel(sprintf('%s',area{1}))
   ylabel(sprintf('%s',area{2}))
   zlabel(sprintf('%s',area{3}))
   set(gca,'CameraPosition',[-11.2424    1.2132    5.1394])
   fp.FormatAxes(gca);  box on; grid on; axis square; 
   fp.FigureSizing(gcf,[3 3 3 3],[10 10 14 10])
      
   %plot in 3D PCA space
   x = cat(3,a{:});
   y = cat(3,b{:});
   x = squeeze(x(:,i,:));   
   y = squeeze(y(:,i,:));
   figure; hold on; 
   x = movmean(x,2,1);
   y = movmean(y,2,1); 
   plot3(x(:,1),x(:,2),x(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',fp.markersizesmall)
   plot3(y(:,1),y(:,2),y(:,3),'marker','o','linestyle','-','linewidth',1,'color',[0.1 0.1 0.8],'MarkerFaceColor',[0.1 0.1 0.8],'MarkerSize',fp.markersizesmall)
   
   %plot a marker at the start
   plot3(x(1,1),x(1,2),x(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   plot3(y(1,1),y(1,2),y(1,3),'marker','d','linestyle','-','linewidth',1,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',fp.markersizesmall)
   xlabel(sprintf('PC %d | %s',i,area{1}))
   ylabel(sprintf('PC %d | %s',i,area{2}))
   zlabel(sprintf('PC %d | %s',i,area{3}))
   set(gca,'CameraPosition',[-18.1734    1.3380    2.5426])
   fp.FormatAxes(gca);  box on; grid on; axis square; 
   fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
   
   
end