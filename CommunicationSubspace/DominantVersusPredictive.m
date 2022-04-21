function [bad_idx,delta]= DominantVersusPredictive(cvl_rrr, cvl_fa,rmv_idx,verbose)
%camden - timesless
%For each predictive dimension in RRR, returns the number of dominant
%dimensions (FAregression) needed to achieve same performance (i.e. within
%SEM across folds). If never reaches, then returns the maximum for models
%tested (usually 10 or max latent Dim of the source population). 

if nargin <3; rmv_idx=[]; end %for plotting a selection of subspaces (i.e. other criteria)
if nargin <4; verbose =1; end

%get the number of dimensions of subspace
rrr_dim = cellfun(@(x) ModelSelect([ mean(x); std(x)/sqrt(size(x,1)) ], 1:size(x,2)), cvl_rrr,'UniformOutput',1);

%get the score of the rrr
score_rrr = cellfun(@(x) 1-x, cvl_rrr,'UniformOutput',0);
score_fa = cellfun(@(x) 1-mean(x), cvl_fa,'UniformOutput',0);

%get the mean - sem for each dimension
min_val = arrayfun(@(n) nanmean(score_rrr{n}(:,1:rrr_dim(n)))-sem(score_rrr{n}(:,1:rrr_dim(n)),1),1:numel(score_rrr),'UniformOutput',0);

dom_d = cell(1,numel(min_val));
for i = 1:numel(min_val)
   temp = arrayfun(@(n) find(score_fa{i}>=min_val{i}(n),1,'first'),1:numel(min_val{i}),'UniformOutput',0);
   dom_d{i} = cat(2,temp{:});
end

%find any times when where the dimensions in rrr_dim > numel(dom_d{i}) and
%replace with min(10,max(fa_dim));
idx = cellfun(@(x,y) numel(x)-numel(y), dom_d,min_val,'UniformOutput',1);
for i = 1:numel(idx)
    if idx(i)==0
    else 
        temp = min(size(cvl_fa{i},2),10);
        dom_d{i} = [dom_d{i},repmat(temp,1,abs(idx(i)))];
    end
end

%get indices for subspaces where any predictive dim is matched with dom
bad_idx = cellfun(@(x) sum(x-(1:size(x,2))==0)>0,dom_d,'UniformOutput',1);

delta = cellfun(@(x) x-(1:size(x,2)),dom_d,'UniformOutput',0);

%plot the figure comparing them
if verbose
   temp_bad_idx = bad_idx(rmv_idx==0); 
   temp_d = cellfun(@(x) [x NaN(1,10-size(x,2))], dom_d,'UniformOutput',0);
   temp_d = cat(1,temp_d{:});  
   temp_d = temp_d(rmv_idx==0,:);
   
   %plot the ones that meet criteria   
   figure; hold on; 
   plot(1:10,1:10,'linestyle','--','linewidth',1.5,'color',[0.5 0.5 0.5])
   xlabel('Number of predictive dimensions')
   ylabel('Minimum number of target dimensions');
   
   a = temp_d(temp_bad_idx==0,:);
   for i = 1:size(a,1)
      y = a(i,~isnan(a(i,:)));
      y = y; %dominant
      x = 1:numel(y); %predictive
      %little jitter for viewing
      y = y+rand(numel(y),1)/4-0.125;
      x = x+rand(numel(x),1)/4-0.125;
      plot(x,y,'marker','o','LineStyle','none','color',[0.5 0.5 0.5],'markersize',3)
   end
   
   %plot the average
   y = nanmean(a);
   y(isnan(y))=[];
   e = sem(a);
   e(isnan(e))=[];
   errorbar(1:numel(y),y,e,'o--','color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',8,'linewidth',1)
   
   
end

end %function
























