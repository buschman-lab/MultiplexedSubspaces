function [rho_all,midx_all] = SubspaceCorr(data_rrr,type)
%Camden MacDowell - timeless
%loop through each target area
verbose=0;
if nargin<2; type = 'spearman'; end
%normalized PEV
n = numel(data_rrr(1).cvl_rrr);
rho_all = cell(1,n);
midx_all = cell(1,n);
for cur_a = 1:n
data = arrayfun(@(n) data_rrr(n).rrr_B{cur_a}, 1:size(data_rrr,2),'UniformOutput',0);

%significance would be by randomly permuting betas within an area I think

%get the rho of for each motif comparison 
ndim = 10;
rho = NaN(size(data,2),size(data,2),ndim);
midx = NaN(size(data,2),size(data,2),ndim);
for i = 1:size(data,2)
    for j = 1:size(data,2)
        if i~=j
           switch type
               case 'spearman'
                   r = corr(data{i}(:,1:ndim),data{j}(:,1:ndim),'type','Spearman'); 
                   [rho(i,j,:),midx(i,j,:)] = max(r,[],2);
               case 'pearson'
                   r = corr(data{i}(:,1:ndim),data{j}(:,1:ndim),'type','Pearson'); 
                   [rho(i,j,:),midx(i,j,:)] = max(r,[],2);                   
               case 'sse'
                   temp = combvec(1:10,1:10)';
                   r = arrayfun(@(n) sse(data{i}(:,temp(n,1)),data{j}(:,temp(n,2))),1:size(temp,1),'UniformOutput',1);
                   r = reshape(r,sqrt(numel(r)),sqrt(numel(r)));
                   [rho(i,j,:),midx(i,j,:)] = min(r,[],2);                   
           end
        else
           rho(i,j,:) = NaN(1,ndim);
           midx(i,j,:) = NaN(1,ndim);
        end
    end    
end

rho_all{cur_a} = rho;
midx_all{cur_a} = midx;

% figure; hold on; 
% t = tiledlayout(4,5);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% for i = 1:ndim
%    nexttile
%    imagesc(rho(:,:,i),[-0.1 0.5]); 
%    axis square; colorbar
%    xlabel('Set dimension of motif'); ylabel('any dimension of motif');
%    title(sprintf('dim %d',i));
% end
% title(t,['Area ',data_rrr(1).area_label{cur_a}],'FontWeight','normal');
% set(gcf,'units','normalized','position',[0 0 1 1])
% 
% %the index of the most similar
% figure; hold on; 
% t = tiledlayout(4,5);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% for i = 1:ndim
%    nexttile
%    imagesc(midx(:,:,i)); 
%    axis square; c=colorbar;
%    ylabel(c,'dimension ID in target');
%    xlabel('set dimension of motif'); ylabel('motif');
%    title(sprintf('dim %d',i));
% end
% title(t,['Area ',data_rrr(1).area_label{cur_a}],'FontWeight','normal');
% set(gcf,'units','normalized','position',[0 0 1 1])



if verbose
    %plot example curve
    t = tiledlayout(5,5);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    count = 1;
    for i = 1:size(data,2)
        if count>25; break; end                
        for j = 1:size(data,2)
            if count>25; break; end                
            if j>i
                nexttile
                hold on;
                r = corr(data{i}(:,1:ndim),data{j}(:,1:ndim)); 
                imagesc(r);
                xlabel(sprintf('Dim Motif %d',i));
                ylabel(sprintf('Dim Motif %d',j));
                axis square
                count = count+1;
                set(gca,'xlim',[0.5 ndim+0.5],'ylim',[0.5 ndim+0.5])
                c=colorbar
                ylabel(c,'rho')
            end        
        end
    end
    set(gcf,'units','normalized','position',[0 0 1 1])
    
end %verbose

end %brain loop region

end %function






