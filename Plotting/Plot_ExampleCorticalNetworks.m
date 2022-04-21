function Plot_ExampleCorticalNetworks(data,cur_rec,motif,area,sigflag,cval)
fp = fig_params_cortdynamics;

data = data{cur_rec};
idx = find(cat(1,data.cur_motif)==motif & cat(1,data.cur_a)==area);

rho_all = data(idx).rho_all; 
sig_thresh = data(idx).sig_thresh; 

%flip rho_all so strongest is positive
for i = 1:size(rho_all,3)
    if nansum(rho_all(:,:,i),'all') < nansum(-1*rho_all(:,:,i),'all')
        rho_all(:,:,i)=-1*rho_all(:,:,i);
    end
end
cone = [0 prctile(rho_all(:),cval)];
%rescale for later dimensions
ctwo = [0 prctile(reshape(rho_all(:,:,2:end),numel(rho_all(:,:,2:end)),1),cval)];

ndim = size(rho_all,3); 
rho_sig = rho_all;
for i = 1:10
   temp = rho_sig(:,:,i);
   temp(abs(temp)<sig_thresh(i))=NaN;
   rho_sig(:,:,i) = temp;
end

%plot networks for each recording | both norm and sig_thresh
for i = 1:ndim
    figure; hold on;
    if i==1
        c = cone;
    else
        c = ctwo;
    end
    if sigflag && sum(~isnan(rho_sig(:,:,i)),'all')>1
        PlotMesoFrame(rho_all(:,:,i)); %get color scaling
%         cval = get(gca,'Clim');
        PlotMesoFrame(rho_sig(:,:,i));
        set(gca,'Clim',c);
        colorbar
        fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    else
        PlotMesoFrame(rho_all(:,:,i));
        set(gca,'Clim',c);
        cc = colorbar;
        set(cc,'ytick',get(cc,'ylim'));
        
        fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    end
    set(gca,'ydir','reverse');
    title(sprintf('motif %d | dim %d | area %d | rec %d',motif, i, area, cur_rec));    
end




end