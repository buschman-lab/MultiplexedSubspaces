function [cv_map,x] = cvStrengthMap(data,type)
%Camden - timeless
%data is a structure with the cv information per motif 
%type is the map that you want to produce

paired_areas = data(1).paired_areas;


switch type
    case 'r_first' %at the maximum timepoint
        z = [data.r];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) x(1),z,'UniformOutput',1);
        
    case 'r_first_norm' %at the maximum timepoint
        z = [data.r];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) x(1)-min(x),z,'UniformOutput',1);        
     
        
    case 'r_all' %average across all significanct and stable timepoints
        z = [data.r];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) nanmean(x),z,'UniformOutput',1);   
        
    case 'r_all_sum' %average across all significanct and stable timepoints
        z = [data.r];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) nansum(x),z,'UniformOutput',1);           
        
    case 'r_first_diag'
        z = [data.r];         
        idx = [data.best_idx];        
        for i = 1:size(idx,1)
            for j = 1:size(idx,2)
                z{i,j}(ismember(idx{i,j},find(diff([data(1).t{1}],1)==0)))=[];
            end
        end   
        
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) x(1),z,'UniformOutput',1); 
        
    case 'r_all_diag' %exclude diagnols        
        z = [data.r];         
        idx = [data.best_idx];        
        for i = 1:size(idx,1)
            for j = 1:size(idx,2)
                z{i,j}(ismember(idx{i,j},find(diff([data(1).t{1}],1)==0)))=[];
            end
        end   
        
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) nanmean(x),z,'UniformOutput',1);      
        
    case 'r_all_diag_sum' %exclude diagnols        
        z = [data.r];         
        idx = [data.best_idx];        
        for i = 1:size(idx,1)
            for j = 1:size(idx,2)
                z{i,j}(ismember(idx{i,j},find(diff([data(1).t{1}],1)==0)))=[];
            end
        end   
        
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) nansum(x),z,'UniformOutput',1);         
        
    case 'aTheta_xv' %average across all significanct and stable timepoints
        z = [data.aTheta_xv];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) (real(x(1))),z,'UniformOutput',1);          
        
    case 'bTheta_xv' %average across all significanct and stable timepoints
        z = [data.bTheta_xv];         
        %fill empty subspaces with NaN
        idx = cellfun(@(x) isempty(x),z);
        z(idx) = {NaN};
        x = cellfun(@(x) (real(x(1))),z,'UniformOutput',1);    
        
    case 'best_idx_one' %average across all significanct and stable timepoints
        z = [data.best_idx];         
        %fill empty subspaces with NaN
        temp_idx = cellfun(@(x) isempty(x),z);
        idx = [data(1).t{1},[NaN;NaN]];
        z(temp_idx) = {size(temp_idx,2)};        
        x = cellfun(@(x) idx(1,x(1)),z,'UniformOutput',1);     
        
    case 'best_idx_two' %average across all significanct and stable timepoints
        z = [data.best_idx];         
        %fill empty subspaces with NaN
        temp_idx = cellfun(@(x) isempty(x),z);
        idx = [data(1).t{1},[NaN;NaN]];
        z(temp_idx) = {size(temp_idx,2)};        
        x = cellfun(@(x) idx(2,x(1)),z,'UniformOutput',1);           
        

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

