function [autorho,arearho,arearho_ordered,areaorder,arearho_avg,inoutrho,inoutrho_avg,xarearho,xarearho_ordered,xarearho_avg] = AnalyzeCorticalNetworks(data,dataout)
%Camden MacDowell - timeless

%autocorrelation within the same subspace (by dimension separation)
autorho = NaN(6,8,14,9);
for cur_rec = 1:6
    for cur_a = 1:8
        for cur_m = 1:14
            motif_list = cat(1,data{cur_rec}.cur_motif);
            area_list = cat(1,data{cur_rec}.cur_a);
            idx = find(motif_list==cur_m & area_list == cur_a);
            if ~isempty(idx)
                x = data{cur_rec}(idx).rho_all;
                x = conditionDffMat(x);
                y = corr(x');
                %do it by dimension (i.e. are similar dimensions/pev similar network?)
                for i = 1:size(x,1)-1    
                    autorho(cur_rec,cur_a,cur_m,i) = nanmean(y((triu(y,i)-triu(y,i+1))~=0)); 
                end
            end
        end
    end
end

%correlation between dimensions across motifs (best match and ordered)
arearho = NaN(6,8,14,10);
areaorder = NaN(6,8,14,10);
arearho_ordered = NaN(6,8,14,10);
arearho_avg = NaN(6,8,14,10);
for cur_rec = 1:6
    for cur_a = 1:8
        motif_list = unique(cat(1,data{cur_rec}.cur_motif));
        for cur_m = 1:numel(motif_list)
            cur_motif = motif_list(cur_m); %contigencies for if there is a missing motif (which doesn't happen in final data but happened when initially writing this. 
            area_list = cat(1,data{cur_rec}.cur_a);
            idx = find(area_list == cur_a);            
            if ~isempty(idx) 
                x = cat(4,data{cur_rec}(idx).rho_all);
                x = arrayfun(@(n) conditionDffMat(x(:,:,:,n)),1:size(x,4),'UniformOutput',0);
                for i = 1:10
                    %best fit                   
                    [rho,ord] = cellfun(@(y) max(abs(corr(x{cur_m}(i,:)',y'))),x,'UniformOutput',1);
                    rho(cur_m)=NaN;
                    ord(cur_m)=NaN;
                    arearho(cur_rec,cur_a,cur_motif,i) = nanmean(rho); 
                    areaorder(cur_rec,cur_a,cur_motif,i) = nanmedian(abs(ord-i)); %get how far away
                    %ordered
                    rho = cellfun(@(y) abs(corr(x{cur_m}(i,:)',y(i,:)')),x,'UniformOutput',1);
                    rho(cur_m)=NaN;
                    arearho_ordered(cur_rec,cur_a,cur_motif,i) = nanmean(rho);     
                    %average
                    rho = cellfun(@(y) nanmean(abs(corr(x{cur_m}(i,:)',y'))),x,'UniformOutput',1);
                    rho(cur_m)=NaN;
                    arearho_avg(cur_rec,cur_a,cur_motif,i) = nanmean(rho);     
                end
            end
        end
    end
end


%correlation between in and out (ordered dimension) | assumes data and
%dataout match (they should) 
%absolute value of the correlation overall sign of the map is abritrary
inoutrho = NaN(6,8,14,10);
inoutrho_avg = NaN(6,8,14,10);
for cur_rec = 1:6
    for cur_a = 1:8
        for cur_m = 1:14
            motif_list = cat(1,data{cur_rec}.cur_motif);
            area_list = cat(1,data{cur_rec}.cur_a);
            idx = find(motif_list==cur_m & area_list == cur_a);
            if ~isempty(idx)
                x = data{cur_rec}(idx).rho_all;
                x = conditionDffMat(x);
                y = dataout{cur_rec}(idx).rho_all;
                y = conditionDffMat(y);
                for i = 1:size(x,1)    
                    inoutrho(cur_rec,cur_a,cur_m,i) = abs(corr(x(i,:)',y(i,:)'));
                    inoutrho_avg(cur_rec,cur_a,cur_m,i) = nanmean(abs(corr(x(i,:)',y'))); %correlation between this motif and all others (to adjust for noise as we get farther along
                end
            end
        end
    end
end

%test correlations between areas (same motif | both ordered and best)
%also needs to be absolute value since the direction is dependent on
%different betas
xarearho = NaN(6,8,14,10);
xarearho_ordered = NaN(6,8,14,10);
xarearho_avg = NaN(6,8,14,10);
for cur_rec = 1:6
    area_list = unique(cat(1,data{cur_rec}.cur_a));
    for cur_a = 1:numel(area_list)
        cur_area = area_list(cur_a);
        for cur_m = 1:14
            motif_list = cat(1,data{cur_rec}.cur_motif);
            idx = find(motif_list==cur_m);
            if ~isempty(idx) 
                x = cat(4,data{cur_rec}(idx).rho_all);
                x = arrayfun(@(n) conditionDffMat(x(:,:,:,n)),1:size(x,4),'UniformOutput',0);
                for i = 1:10
                    %best fit
                    rho = cellfun(@(y) max(abs(corr(x{cur_a}(i,:)',y'))),x,'UniformOutput',1);
                    rho(cur_a)=NaN;
                    xarearho(cur_rec,cur_area,cur_m,i) = nanmean(rho); 
                    %ordered
                    rho = cellfun(@(y) abs(corr(x{cur_a}(i,:)',y(i,:)')),x,'UniformOutput',1);
                    rho(cur_a)=NaN;
                    xarearho_ordered(cur_rec,cur_area,cur_m,i) = nanmean(rho);  
                    rho = cellfun(@(y) nanmean(abs(corr(x{cur_a}(i,:)',y'))),x,'UniformOutput',1);
                    rho(cur_a)=NaN;
                    xarearho_avg(cur_rec,cur_a,cur_motif,i) = nanmean(rho);                        
                end
            end
        end
    end
end


end %function 































