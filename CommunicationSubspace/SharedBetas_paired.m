function [rsq,matchdim] = SharedBetas_paired(data,type,ndim,shufseed,matchdim)
if nargin <5; matchdim=[]; end
%Camden - timeless

area_label = data(1).area_label;
%loop through each area
if isempty(matchdim)
    matchdim = NaN(numel(area_label),size(data,2),size(data,2),ndim);
    rsq = NaN(numel(area_label),size(data,2),size(data,2),ndim);
else
    rsq = NaN(numel(area_label),size(data,2),size(data,2));
end
for cur_area = 1:numel(area_label)
    switch type
%for each pairing with beta weights in area X. Use the 
        case 'projection'
            fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));
            for cur_motif = 1:size(data,2)
                %get the activity of the target area
                x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});
                %normalize to baseline
                x = normalizeToBaseline(x,[1:2],'mean');
                %use post stimulus
                x = x(:,3:end,:);
                %remove the psth
                x = x-nanmean(x,3);
                x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
                for comp_motif = 1:size(data,2)        
                    Y = cell(1,ndim);
                    Yhat = cell(1,ndim);
                    if comp_motif == cur_motif
                    else
                        for cur_d = 1:ndim
                            %get the predicted activity for the target fit
                            B = data(cur_motif).rrr_B{cur_area}(:,cur_d);
                            V = data(cur_motif).rrr_V{cur_area}(:,cur_d);
                            Y{cur_d} = x*B*V';

                            %get the best dim in target subspace
                            BB = data(comp_motif).rrr_B{cur_area}(:,cur_d);
                            if shufseed>1
                                rng(shufseed)
                                shufidx = ShuffleWithinArea(BB,data,cur_motif,cur_area);
                            else
                                shufidx = repmat(1:size(BB,1),size(BB,2),1)';    
                            end
                            VV = data(comp_motif).rrr_V{cur_area}(:,cur_d);
                            Yhat{cur_d} = x*BB(shufidx)*VV';
                        end
                        idx = combvec(1:ndim,1:ndim);
                        rtemp = arrayfun(@(n) 1-NormalizedSquaredError(Y{idx(2,n)}, Yhat{idx(1,n)}), 1:size(idx,2),'UniformOutput',1);
                        %for each dimension in Y, get the best from Yhat
                        [rsq(cur_area,cur_motif,comp_motif,1:ndim),tempmatchdim] = arrayfun(@(n) max(rtemp(idx(2,:)==n)), unique(idx(2,:)),'UniformOutput',1);
                        matchdim(cur_area,cur_motif,comp_motif,1:ndim) = tempmatchdim;
                    end %if else
                end %comp motif
            end %motif
        case 'projection_reverse'
            fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));
            for cur_motif = 1:size(data,2)
                %get the activity of the target area
                x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==1});
                %normalize to baseline
                x = normalizeToBaseline(x,[1:2],'mean');
                %use post stimulus
                x = x(:,3:end,:);
                %remove the psth
                x = x-nanmean(x,3);
                x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
                for comp_motif = 1:size(data,2)        
                    Y = cell(1,ndim);
                    Yhat = cell(1,ndim);
                    if comp_motif == cur_motif
                    else
                        for cur_d = 1:ndim
                            %get the predicted activity for the target fit
                            B = data(cur_motif).rrr_B{cur_area}(:,cur_d);
                            V = data(cur_motif).rrr_V{cur_area}(:,cur_d);
                            Y{cur_d} = x*B*V';

                            %get the best dim in target subspace
                            BB = data(comp_motif).rrr_B{cur_area}(:,cur_d);
                            if shufseed>1
                                rng(shufseed)
                                shufidx = ShuffleWithinArea(BB,data,cur_motif,cur_area);
                            else
                                shufidx = repmat(1:size(BB,1),size(BB,2),1)';    
                            end
                            VV = data(comp_motif).rrr_V{cur_area}(:,cur_d);
                            Yhat{cur_d} = x*BB(shufidx)*VV';
                        end
                        idx = combvec(1:ndim,1:ndim);
                        rtemp = arrayfun(@(n) 1-NormalizedSquaredError(Y{idx(2,n)}, Yhat{idx(1,n)}), 1:size(idx,2),'UniformOutput',1);
                        %for each dimension in Y, get the best from Yhat
                        [rsq(cur_area,cur_motif,comp_motif,1:ndim),tempmatchdim] = arrayfun(@(n) max(rtemp(idx(2,:)==n)), unique(idx(2,:)),'UniformOutput',1);
                        matchdim(cur_area,cur_motif,comp_motif,1:ndim) = tempmatchdim;
                    end %if else
                end %comp motif
            end %motif
        case 'fullprojection'            
            fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));
            for cur_motif = 1:size(data,2)
                %get the activity of the target area
                x = cat(1,data(cur_motif).area_val{ismember(1:numel(area_label),cur_area)==0});
                %normalize to baseline
                x = normalizeToBaseline(x,[1:2],'mean');
                %use post stimulus
                x = x(:,3:end,:);
                %remove the psth
                x = x-nanmean(x,3);
                x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
                for comp_motif = 1:size(data,2)        
                    if comp_motif == cur_motif
                    else
                        tempmatchdim = matchdim(cur_area,cur_motif,comp_motif,1:ndim);
                        %get the predicted activity for the target fit
                        B = data(cur_motif).rrr_B{cur_area}(:,1:ndim);
                        V = data(cur_motif).rrr_V{cur_area}(:,1:ndim);
                        Y = x*B*V';

                        %get the best dim in target subspace
                        BB = data(comp_motif).rrr_B{cur_area}(:,tempmatchdim);
                        if shufseed>1
                            rng(shufseed)
                            shufidx = ShuffleWithinArea(BB,data,cur_motif,cur_area);
                        else
                            shufidx = repmat(1:size(BB,1),size(BB,2),1)';    
                        end
                        VV = data(comp_motif).rrr_V{cur_area}(:,tempmatchdim);
                        Yhat= x*BB(shufidx)*VV';
                        rsq(cur_area,cur_motif,comp_motif) = 1-NormalizedSquaredError(Y,Yhat);
                    end %if else
                end %comp motif
            end %motif
    otherwise
        error('unknown type')
    end
end %area




end %function 


function shufidx = ShuffleWithinArea(B,data,cur_motif,cur_area)

grp = data(cur_motif).grouping{cur_area};
if numel(grp)~=size(B,1) %working on the single area
   grp = ones(size(B,1),1);
end
ugrp = unique(grp);

shufidx = NaN(size(B));
for j = 1:size(B,2)
    tempidx = 1:numel(grp);
    for i = 1:numel(ugrp)
        temp = tempidx(grp==ugrp(i));
        %maintain sign of neuron
        s = B(grp==ugrp(i),j)>0;
        a = temp(s==1);
        temp(s==1) = a(randperm(numel(a),numel(a)));
        a = temp(s==0);
        temp(s==0) = a(randperm(numel(a),numel(a)));
        tempidx(grp==ugrp(i)) = temp;
    end
    shufidx(:,j) = tempidx;
end

end






















