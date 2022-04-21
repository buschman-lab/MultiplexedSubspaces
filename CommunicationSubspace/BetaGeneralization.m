function [rsq_full,rsq_byarea] = BetaGeneralization(data_rrr,match_best)
%Camden - timeless

if nargin <2; match_best=1; end

if match_best==1
   [~,midx_all] = SubspaceCorr(data_rrr,'sse'); 
end
area_val = data_rrr(1).area_val;
area_label = data_rrr(1).area_label;
paired_areas = data_rrr(1).paired_areas;

rsq_full = cell(1,numel(area_label));
rsq_byarea = cell(1,numel(area_label));
for i = 1:numel(area_label)
    grp = data_rrr(1).grouping{i};
    idx = strcmp(area_label,area_label{paired_areas(i)});
    x = cat(1,area_val{idx==0,:});    

    %normalize to baseline
    x = normalizeToBaseline(x,[1:2],'mean');

    %use post stimulus
    x = x(:,3:end,:);

    %get the psth
    x = x-nanmean(x,3);
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    
    if match_best
        best_dim = midx_all{i};
    else
        best_dim = arrayfun(@(n) n*ones(size(data_rrr,2)),1:10,'UniformOutput',0);
        best_dim = cat(3,best_dim{:});
    end
    
    %full model loop through each pair of motifs across all dimensions   
    for ii = 1:size(data_rrr,2)
        for jj = 1:size(data_rrr,2)
            if ii~=jj
                B = data_rrr(ii).rrr_B{i}(:,1:size(best_dim,3));
                V = data_rrr(ii).rrr_V{i}(:,1:size(best_dim,3)); 
                Y = x*B*V';

                %get the best dim in target subspace
                BB = data_rrr(jj).rrr_B{i}(:,best_dim(ii,jj,:));
                VV = data_rrr(jj).rrr_V{i}(:,best_dim(ii,jj,:));
                Yhat = x*BB*VV';

                rsq_full{i}(ii,jj) = 1-NormalizedSquaredError(Y, Yhat);
            else
                rsq_full{i}(ii,jj) = NaN;
            end
        end
    end    
    
    
%     %full model loop through each pair of motifs    
%     for ii = 1:size(data_rrr,2)
%         for jj = 1:size(data_rrr,2)
%             if ii~=jj
%                 for d = 1:size(best_dim,3)
%                     B = data_rrr(ii).rrr_B{i}(:,d);
%                     V = data_rrr(ii).rrr_V{i}(:,d); 
%                     Y = x*B*V';
% 
%                     %get the best dim in target subspace
%                     BB = data_rrr(jj).rrr_B{i}(:,best_dim(ii,jj,d));
%                     VV = data_rrr(jj).rrr_V{i}(:,best_dim(ii,jj,d));
%                     Yhat = x*BB*VV';
% 
%                     rsq_full{i}(ii,jj,d) = 1-NormalizedSquaredError(Y, Yhat);
%                 end
%             else
%                 rsq_full{i}(ii,jj,1:size(best_dim,3)) = NaN(1,size(best_dim,3));
%             end
%         end
%     end    
    
    %loop through each brain region
    unique_grp = unique(grp);
    temp = NaN(1,numel(area_label));
    for cur_grp = 1:numel(unique_grp)        
        idx = grp==unique_grp(cur_grp);
        Y = x(:,idx)*B(idx,:)*V';
        Yhat = x(:,idx)*BB(idx,:)*VV';
        temp(1,unique_grp(cur_grp)) = 1-NormalizedSquaredError(Y, Yhat);        
    end
    rsq_byarea{i} = temp;

end %subspace identification loop

end %function 

























