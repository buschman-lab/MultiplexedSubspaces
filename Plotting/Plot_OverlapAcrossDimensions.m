function Plot_OverlapAcrossDimensions(data_all,method)
%Computes the spatial overlap of subspace networks within a subspace

%%
if nargin < 2; method = 1; end %def(1) use correlation 2=use overlap

% %% get statistics for specific overlaps used in text
% %rec 1, motif 8, area 2 (motor)
% temp = ComputeOverLap(data_all{1},areas(2),8,method,1);
% temp = reshape(temp,10,10);
% 
% %get the overlap between dimensions 4 and 8
% a = temp(8,4)*100;
% b = temp(9,6)*100;
% c = nanmean([temp(9,8),temp(9,4),temp(8,6),temp(6,4)])*100;
% 
% rho = numPxls(data_all{1},2,8);
% x = rho(:,:,[4,8,6,9]);
% n = sum(~isnan(reshape(x,size(x,1)*size(x,2),size(x,3))));
% n = n(1);
% x = nansum(reshape(x,size(x,1)*size(x,2),size(x,3)))/n;
% %%
% y = binocdf(a*n,n,(x(1)*x(2)),'upper');
% y = binocdf(b*n,n,(x(3)*x(4)),'upper');
% %between
% 1-binocdf(temp(6,4)*n,n,(x(1)*x(3)),'upper');
% 1-binocdf(temp(9,4)*n,n,(x(1)*x(4)),'upper');
% 1-binocdf(temp(8,6)*n,n,(x(2)*x(3)),'upper');
% 1-binocdf(temp(9,8)*n,n,(x(2)*x(4)),'upper');

%%
ovrlap = NaN(6,8,14,45);
for cur_rec = 1:6
    data = data_all{cur_rec}; 
    areas = unique(cat(1,data.cur_a));
    for cur_area = 1:numel(areas)
        for cur_motif = 1:14
           ovrlap(cur_rec,areas(cur_area),cur_motif,:) = ComputeOverLap(data,areas(cur_area),cur_motif,method);
        end
    end
end

%% do full histogram and overlay lines for each brain area. 
fp = fig_params_cortdynamics;

if method ==1
    binsize = 0.1;
    edges = -1:binsize:1;
else
    binsize = 0.1;
    edges = 0:binsize:1;
end
N = NaN(numel(edges)-1,8);
for i = 1:8
    x = squeeze(ovrlap(:,i,:,:));
    N(:,i) = histcounts(x(:),edges);
    N(:,i) = N(:,i)/sum(N(:,i));
end

%full
yall = histcounts(ovrlap(:),edges);
yall = yall/sum(yall);
yall = smoothdata(yall,'gaussian',5);

%fraction with <10% overlap
numel(find(ovrlap(:)<.1))/numel(ovrlap(:))

figure; hold on; 
histogram(ovrlap(:),edges,'Normalization','probability','FaceAlpha',0.25,'FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.25)
col = fp.c_area;
for i = 1:8
    if method ==1
        y = smoothdata(N(:,i),'gaussian',10);
    else
%         y = smoothdata(N(:,i),'gaussian',3);
        y = N(:,i);
    end
    plot(edges(1:end-1)+binsize/2,y,'linewidth',1,'color',col(i,:),'linestyle','-')
end
% plot(edges(1:end-1)+binsize/2,yall,'linewidth',3,'color',[0.5 0.5 0.5])

ylabel('Probability');
xlabel('Subspace Network Similarity');
fp.FormatAxes(gca); box on; grid off
fp.FigureSizing(gcf,[3 2 3 3],[2 10 10 10])

ymax = get(gca,'ylim');
xavg = nanmean(ovrlap(:));
plot([xavg,xavg],[0,ymax(2)],'linewidth',1,'color','r','linestyle','--');
plot([0,0],[0,ymax(2)],'linewidth',1,'color','k','linestyle','--');

% ci = bootci(1000,@nanmean,ovrlap(:))
title(sprintf('mu %0.4f CI = %0.4f',xavg,std(ovrlap(:))),'fontweight','normal');

end


function ovrlap = ComputeOverLap(data,cur_area,cur_motif,method,flat)
    if nargin <5; flat=0; end
    
    %get the motif and area index
    idx = find(cat(1,data.cur_motif)==cur_motif & cat(1,data.cur_a)==cur_area);

    rho_all = data(idx).rho_all; 
    %direction of correlation is arbitrary. Flip so strongest are positive -
    %that is our subspace network for this dimension. 
    for i = 1:size(rho_all,3)
        if nansum(rho_all(:,:,i),'all') < nansum(-1*rho_all(:,:,i),'all')
            rho_all(:,:,i)=-1*rho_all(:,:,i);
        end
    end
    if method == 1
        rho = rho_all; 
    else
        %parse data structure
        rho = rho_all; 
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
    end
    
    %get overlap 
    ovrlap = NaN(size(rho,3),size(rho,3));
    for i = 1:size(rho,3)
        for j = 1:size(rho,3)
            if i>j %upper triangle
                a = squeeze(rho(:,:,i));
                b = squeeze(rho(:,:,j));
                %overlap is 2x the number of shared pixels over the sum of
                %the two subnetowrks total pixels
                if method==1
                    ovrlap(i,j) = corr(a(:),b(:),'rows','complete');
                else
                    c = a+b;
                    ovrlap(i,j) = (2*sum(c==2))/(sum(a==1)+sum(b==1));
                end                
            end            
        end
    end      
    
    if flat
        ovrlap = ovrlap(:);
    else
        ovrlap = ovrlap(~isnan(ovrlap));
    end
    
end


function rho = numPxls(data,cur_area,cur_motif)
  
    %get the motif and area index
    idx = find(cat(1,data.cur_motif)==cur_motif & cat(1,data.cur_a)==cur_area);

    rho_all = data(idx).rho_all; 
    %direction of correlation is arbitrary. Flip so strongest are positive -
    %that is our subspace network for this dimension. 
    for i = 1:size(rho_all,3)
        if nansum(rho_all(:,:,i),'all') < nansum(-1*rho_all(:,:,i),'all')
            rho_all(:,:,i)=-1*rho_all(:,:,i);
        end
    end
    %parse data structure
    rho = rho_all; 
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

    
end





