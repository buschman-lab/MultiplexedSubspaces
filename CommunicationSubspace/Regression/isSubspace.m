function idx = isSubspace(data, verbose)
% Camden - timeless
% returns the indices of subspaces that meet the criteria to be a subspace
% @Criteria based off Semedo et al., Neuron 2019 | Cortical Areas Interact
% through a communication subspace
% @inputs: data is a structure with cvl_fa, cvl_ridge, cv_rrr, etc. from
% RRR. If data contains multiple rows, will loop through each


%get the appropriate full model 
[score,idx] = cellfun(@(x) bestLambda(x), data(cur_fit).cvl_ridge,'UniformOutput',0);

%get the number of dimensions of subspace
rrr_dim = cellfun(@(x) ModelSelect([ mean(x); std(x)/sqrt(size(x,1)) ], 1:size(x,2)), data(cur_fit).cvl_rrr,'UniformOutput',1);

%get the score of the first rrr dimension 
score_rrr = cellfun(@(x) 1-x(:,1), data(cur_fit).cvl_rrr,'UniformOutput',0);

%error from full model relative to the total performance
[e_rat,~,~,withinSEM] = cellfun(@(x,y) errBetween(x,y), score,data(cur_fit).cvl_rrr,'UniformOutput',1);

%Criteria 1: full model is nonpredictive (xval produced negative) or fully predictive (self)
y = cellfun(@nanmin,score);
bad_idx = y<0 | y==1;

%Criteria 2: first dimension of RR has to be predictive (i.e. all cross
%validations above zero performance)
y = cellfun(@nanmin,score_rrr);
bad_idx = bad_idx | y<=0; 

%Criteria 3: rrr_dim must be smaller dimensionality of shared variance in the local population
y = cat(1,data(cur_fit).qOpt{:}) - rrr_dim;
bad_idx = bad_idx | y<=0;
 
%Criteria 4: rrr_dim must be smaller dimensionality of shared variance in the target population
% y = cat(1,data(cur_fit).qOpt_target{:}) - rrr_dim;
% bad_idx = bad_idx | y<=0;

%parse the 'strength' by how much of full model is explained
strength = e_rat; strength(e_rat<=0.1)=1; 
strength(e_rat<=0.25 & e_rat>=0.1)=2; 
strength(e_rat<=1 & e_rat>=0.25)=3;

if verbose
   close all
   arrayfun(@(n) plot_rrrSummary(score{n},data(cur_fit).cvl_rrr{n},data(cur_fit).cvl_fa{n},bad_idx(n),strength(n)),1:numel(idx))
   
   %plot the distribution of model strength
   e_rat(bad_idx)=[]; withinSEM(bad_idx)=[]; 
   figure; hold on; 
   histogram(e_rat(withinSEM==1),'binwidth',0.01,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1]);
   histogram(e_rat(withinSEM==0),'binwidth',0.01,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8]);
   legend('Within SEM','Outside SEM','location','northeastoutside'); xlabel('performance fraction difference'); ylabel('# of subspaces')
   
   %confirm that 'strong' is not just really weakly predictive ones
   temp = score; 
   temp(bad_idx)=[];
   temp = cellfun(@nanmean,temp); 
   figure; hold on;   
   histogram(temp(withinSEM==1),'binwidth',0.05,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1]);
   histogram(temp(withinSEM==0),'binwidth',0.05,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8]);
   legend('Within SEM','Outside SEM','location','northeastoutside'); xlabel('performance fraction difference'); ylabel('# of subspaces')      
end



% for cur_fit = 1:size(data,1); 



end %function end





