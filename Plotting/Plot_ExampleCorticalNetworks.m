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
        fH = PlotMesoFrame(rho_all(:,:,i)); %get color scaling
        
%         x = ~isnan(rho_sig(:,:,i));
%         x = imfill(x,8,'holes');
%         b = bwboundaries(x);
%         
%         for k=1:numel(b)
%            plot(b{k}(:,2),b{k}(:,1),'color',fp.c_lr,'LineWidth',1);
%         end
        %change the alpha
        
        for k = 1:2
            x = ~isnan(rho_sig(:,:,i));
            y = get(fH{k},'AlphaData');
            mask = y; 
            mask(y==1)=0.60;
            mask(x==1 & y==1)=1;        
            set(fH{k},'AlphaData',mask)        
        end     
        
%         cval = get(gca,'Clim');
        %Plot an outline
%         PlotMesoFrame(rho_sig(:,:,i));
        set(gca,'Clim',c);
        colorbar
        fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
        cc = colorbar;
        set(cc,'ytick',get(cc,'ylim'));
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