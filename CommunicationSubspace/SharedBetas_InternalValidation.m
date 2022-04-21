function rsq = SharedBetas_InternalValidation(data,reverseFlag)
%Camden - timeless
if nargin <2; reverseFlag=0; end
ndim = 10;
nxval = 10;
rng('default')
area_label = data(1).area_label;
%loop through each area
rsq = NaN(numel(area_label),size(data,2),ndim);
for cur_area = 1:numel(area_label) 
    fprintf('\n\t working on area %d of %d',cur_area,numel(area_label));
    for cur_motif = 1:size(data,2)
        %split trials in half
        %get the activity of the target area
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
        x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';             
        y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';        
        
        %10 fold cross validation        
        cvp = cvpartition(size(x,1),'kfold',nxval);
        
        %get the lambda (using full rec as done with initial fitting)
        dMaxShrink = .5:.01:1;
        lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);        
        [~,idx] = bestLambda(data(cur_motif).cvl_ridge{cur_area});
        
        % Cross validated Betas
        B = cell(1,nxval); V = cell(1,nxval);
        for i = 1:nxval
            [~,B{i},V{i}] = ReducedRankRegress(y(cvp.training(i),:), x(cvp.training(i),:), ndim,'scale',false,'RIDGEINIT',lambda(idx));               
        end
        % Full betas
        [~,BB,VV] = ReducedRankRegress(y, x, ndim,'scale',false,'RIDGEINIT',lambda(idx));  
        
        %get the original index
        fullidx = arrayfun(@(n) find(cvp.test(n)==1),1:nxval,'UniformOutput',0);
        fullidx = cat(1,fullidx{:});
        for cur_d = 1:ndim
            %project with the withheld
            Ytemp = arrayfun(@(n) x(cvp.test(n),:)*B{n}(:,cur_d)*V{n}(:,cur_d)',1:nxval,'UniformOutput',0);
            Yhat = y;
            Yhat(fullidx,:) = cat(1,Ytemp{:});
            
            %Get the true full
            Y = x*BB(:,cur_d)*VV(:,cur_d)';
            rsq(cur_area,cur_motif,cur_d) = 1-NormalizedSquaredError(Y,Yhat);
        end
    end
end %area


end %function 

%reason for doing this, instead of split haves is the number of
%observations x predictors... basically, this is more fair

% x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';     
% xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';  
% y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';    
% yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';    
% x_full = reshape(x_full,[size(x_full,1),size(x_full,2)*size(x_full,3)])';





















