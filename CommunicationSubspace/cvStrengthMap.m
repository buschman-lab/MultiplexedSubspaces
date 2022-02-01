function cv_map = cvStrengthMap(data,motif_cvs,type)
%Camden - timeless
%data is a structure with the cv information per motif 
%type is the map that you want to produce

paired_areas = data(1).paired_areas;


switch type
    case 'r'
        z = [data.r];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        
        %do the same for the indices of each motif
        y = motif_cvs;
        
        x = NaN(size(y)); 
        for i = 1:size(y,1)
            for j = 1:size(y,2)
               if ~isempty(y{i,j})
               x(i,j) = sum(z{i}(y{i,j}));
               end
            end
        end
        
    case 'r_norm'
        z = [data.r_norm];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        
        %do the same for the indices of each motif
        y = motif_cvs;
        
        x = NaN(size(y)); 
        for i = 1:size(y,1)
            for j = 1:size(y,2)
               if ~isempty(y{i,j})
               x(i,j) = sum(z{i}(y{i,j}));
               end
            end
        end      
        
%     case 'pval'        
%         x = [data.pval];         
%         %fill empty subspaces with NaN
%         idx = cellfun(@(x) isempty(x),x);
%         x(idx) = {NaN};
%         x = cellfun(@(x) x(:,1),x); %strength of the first CV  
        
%     case 'cv_weight_norm'
%         y = [data.U];  
%         yy = [data.V];
%         x = NaN(size(y));
%         for idx = 1:size(y,1)                  
%             for j = 1:size(y,2)
%                if ~isempty(y{idx,j})
%                    temp = y{idx,j}(:,1).*yy{idx,j}(:,1);
%                    temp = reshape(temp,[21,size(temp,1)/21]);
%                    temp = zscore(temp);
%                    x(idx,j) = max(nanmean(temp,2));
%                end
%             end
%         end
            
        
%     case 'pca_theta'
%         %for each subspace, get the average angle in both regions to pca
%         x = squeeze(nanmean(cat(3,data.pca_theta),2));     
    case 'num_neu' %a gut check
        x = NaN(size(paired_areas,1),1);
        for i = 1:size(paired_areas,1)
            x(i) = size(data(1).area_val{paired_areas(i,1)},1) + size(data(1).area_val{paired_areas(i,2)},1);
        end
                        
    otherwise
        error('unknown type');
end

%%
cv_map = statmap(x,paired_areas);

end %function end

function y = statmap(x,paired_areas)
   n = max(unique(paired_areas));
   y = NaN(n,n,size(x,2));
   for cur_m = 1:size(x,2)
       for i = 1:size(paired_areas,1)
           y(paired_areas(i,1),paired_areas(i,2),cur_m)= x(i,cur_m);
       end
   end
end

