function [rho,dom,pev,ex,rel_rsq_all] = SubspaceDim_Trial(data,reverseFlag)
%Camden - timeless
%for each trial project activity along a dimension and see how much of the
%full PEV is captured by that dimension. 

if nargin <2; reverseFlag=0; end
ndim = 15;
rng('default')
area_label = data(1).area_label;

%loop through each area
rho = NaN(numel(area_label),size(data,2));
dom = NaN(numel(area_label),size(data,2),ndim);
pev = NaN(numel(area_label),size(data,2),ndim);
rel_rsq_all = cell(numel(area_label),size(data,2));
for cur_area = 1:numel(area_label) 
    fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));    
    for cur_motif = 1:size(data,2)
        %split trials in half
        if reverseFlag ==1
            x = data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==1};
            y = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});            
        else
            x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});
            y = data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==1};
        end
        %normalize to baseline
        x = normalizeToBaseline(x,[1:2],'mean');
        y = normalizeToBaseline(y,[1:2],'mean');
        %use post stimulus
        x = x(:,3:end,:);
        y = y(:,3:end,:);
        %remove the psth
        x = x-nanmean(x,3);
        y = y-nanmean(y,3);          

        %full model
        B = data(cur_motif).rrr_B{cur_area};
        V = data(cur_motif).rrr_V{cur_area};        
        
        rsq = NaN(size(x,3),ndim);
        p = NaN(1,ndim);
        for t = 1:size(x,3)
            xx = squeeze(x(:,:,t))';             
            
            %total predicted activity across all dimensions
            yy = xx*B*V';
            
            %explainable variance of that trial
            p(t) = 100*(1-NormalizedSquaredError(squeeze(y(:,:,t))',yy));

            %how much of PEV captured by each dim PER Trial
            for cur_d = 1:ndim
                Yhat = xx*B(:,1:cur_d)*V(:,1:cur_d)'; %Note: this should be cummulative sum (i.e., same as initial testing)
                rsq(t,cur_d) = 1-NormalizedSquaredError(yy,Yhat);
            end
        end
        %set negatives to zero (so don't inflate subsequent values) and get
        %relative contribution
        rsq(rsq<0)=0; 
        rel_rsq = cat(2,rsq(:,1), diff(rsq,[],2));
        rel_rsq_all{cur_area,cur_motif} = rel_rsq;
        
        %average spearman correlation
        r = triu(corr(rel_rsq','type','spearman'),1);
        r(r==0)=NaN;
        rho(cur_area,cur_motif) = nanmean(r(:));
        
        %dominating dimensions for each trial
        [~,d] = max(rel_rsq,[],2);
        dtemp = arrayfun(@(n) sum(d==n),unique(d),'UniformOutput',1);
        dom(cur_area,cur_motif,1:numel(dtemp))=dtemp;
        
        %explainable variance for trials dominated by each dimension
        ptemp = arrayfun(@(n) nanmean(p(d==n)),unique(d),'UniformOutput',1);
        pev(cur_area,cur_motif,1:numel(ptemp))=ptemp;
        
        %save off all infromation for supplemental figure
        if cur_area == 8 & cur_motif==5
            ex.rsq = rel_rsq;
            ex.rho = r;
            ex.d = d;
            ex.pev = pev;
        end  
    end %motif loop

end %area





end %function 





%% yeah this is not appropriate because they SHOULD be correlated
%     %average spearman correlation across trials in different motifs
%     for cur_motif = 1:size(data,2)
%         temp_rho = NaN(1,size(data,2));
%         for targ_motif = 1:size(data,2)
%             if cur_motif~=targ_motif
%                 r = triu(corr(rel_rsq_all{cur_motif}',rel_rsq_all{targ_motif}','type','spearman'),1);
%                 r(r==0)=NaN;
%                 temp_rho(targ_motif) = nanmean(r(:));   
%                 if cur_area==8 & cur_motif==5
%                     ex.rnull{targ_motif} = r(:);
%                 end
%             end
%         end
%         rho_null(cur_area,cur_motif) = nanmean(temp_rho);
%     end








