function Figure6a(data_all)
%Computes the spatial overlap of subspace networks within a subspace

ovrlap = NaN(6,8,14,68*68);
for cur_rec = 1:6
    data = data_all{cur_rec}; 
    areas = unique(cat(1,data.cur_a));
    for cur_area = 1:numel(areas)
        for cur_motif = 1:14
           ovrlap(cur_rec,areas(cur_area),cur_motif,:) = ComputeNumDimParticipate(data,areas(cur_area),cur_motif);
        end
    end
end

%% do full histogram and overlay lines for each brain area. 
fp = fig_params_cortdynamics;

y = ovrlap; 
y = reshape(y,size(y,1)*size(y,2)*size(y,3),size(y,4));
y(sum(isnan(y),2)==size(y,2),:)=[];
%keep mask of all animals
y(:,sum(isnan(y),1)>0)=NaN;
y = nanmean(y);
y = reshape(y,[68,68]);
%make into fraction
y = y/10;


PlotMesoFrame(y);
set(gca,'Clim',[0,.4]);
cc = colorbar;
set(cc,'ytick',get(cc,'ylim'));
fp.FigureSizing(gcf,[3 2 6 6],[10 10 12 11])
title('Network Distribution','fontweight','normal');

end


function y = ComputeNumDimParticipate(data,cur_area,cur_motif)
    %get the motif and area index
    idx = find(cat(1,data.cur_motif)==cur_motif & cat(1,data.cur_a)==cur_area);

    rho = data(idx).rho_all; 
    sig_thresh = data(idx).sig_thresh; 
    %threshold and binarize
    for i = 1:size(rho,3)
       temp = rho(:,:,i);
       temp(abs(temp)<sig_thresh(i))=0;
       rho(:,:,i) = temp;
    end
    %binarize. anything that is significantly positively correlated is our network
    rho(rho>0)=1; 
    rho(rho<=0)=0;      
    
    y = sum(rho,3);
    y = y(:);
end






