function [good_idx,stats,bad_idx_all, label] = isSubspace(data, verbose)
% Camden - timeless
% returns the indices of subspaces that meet the criteria to be a subspace
% @Criteria based off Semedo et al., Neuron 2019 | Cortical Areas Interact
% through a communication subspace
% @inputs: data is a structure with cvl_fa, cvl_ridge, cv_rrr, etc. from
% RRR. If data contains multiple rows, will loop through each

good_idx = cell(1,size(data,2));
label=[];
stats = struct();
bad_idx_all = cell(1,size(data,2));
for cur_fit = 1:size(data,2)
    %get the appropriate full model 
    score = cellfun(@(x) bestLambda(x), data(cur_fit).cvl_ridge,'UniformOutput',0);

    %get the number of dimensions of subspace
    rrr_dim = cellfun(@(x) ModelSelect([ mean(x); std(x)/sqrt(size(x,1)) ], 1:size(x,2)), data(cur_fit).cvl_rrr,'UniformOutput',1);

    %get the score of the rrr
    score_rrr = cellfun(@(x) 1-x, data(cur_fit).cvl_rrr,'UniformOutput',0);

    %error from full model relative to the total performance
    [e_rat,~,~,withinSEM,mean_score] = cellfun(@(x,y) errBetween(x,y), score,data(cur_fit).cvl_rrr,'UniformOutput',1);

    %Criteria 1: full model is nonpredictive (xval produced negative) or fully predictive (self)
    y = cellfun(@nanmin,score);    
    label = {'self'};        
    bad_idx = y==1;
    label(end+1) = {'negative full model'}; 
    bad_idx = [bad_idx, y<0];  %nonpredictive
    

    %Criteria 2: dimensions of RR have to be predictive (i.e. above zero performance)
    y = cellfun(@(x) nanmin(nanmean(x)),score_rrr);
    label(end+1) = {'rrr_predictive'};    
    bad_idx = [bad_idx, y<=0]; 
    
    %Criteria 3: the dominant dimensions on source must be different 
    %and greater in number than the predictive subspace dimensions
    temp_bad_idx = DominantVersusPredictive(data(cur_fit).cvl_rrr, data(cur_fit).cvl_fa,[],0);
    label(end+1) = {'private_dim'};    
    bad_idx = [bad_idx, temp_bad_idx'];

    %Criteria 4: rrr_dim must be smaller dimensionality of shared variance in the target population
    y = cat(1,data(cur_fit).qOpt_target{:}) - rrr_dim;
    label(end+1) = {'target_dim'};
    bad_idx = [bad_idx, y<=0];
    
    %Criteria 5: betas are not driven by single neurons (i.e. anatomical
    %bottleneck). Look at the first dim because this is the one driving the most 
    %variance (i.e. if it passes mustard, then it's okay to have subsequent
    %that are donimanted by one, since the subspace is the combination of all.     
    y = cellfun(@(x) sum(abs(x)>((prctile(abs(x),75)+iqr(abs(x))*5)))>0, data(cur_fit).rrr_B,'UniformOutput',0);
%     y = arrayfun(@(n) sum(y{n}(1:rrr_dim(n))), 1:numel(y))';
    y = arrayfun(@(n) sum(y{n}(1)), 1:numel(y))';
    label(end+1) = {'beta_outlier'};
    bad_idx = [bad_idx, y>0];
    
    %Criteria 6: optimal RRR has to reach at least 0.70 of the full model
    %variance
    label(end+1) = {'poor_fit'};
    bad_idx = [bad_idx, e_rat<0.7];
    
    %compile into a single good index
    good_idx{cur_fit} = sum(bad_idx,2)==0;
    bad_idx_all{cur_fit} = bad_idx;
    bad_idx = sum(bad_idx,2)>0;

    stats(cur_fit).full_score = score;
    stats(cur_fit).rrr_score = arrayfun(@(n) score_rrr{n}(:,rrr_dim(n)),1:numel(rrr_dim),'UniformOutput',0)';
    stats(cur_fit).ss_dim = rrr_dim;
    stats(cur_fit).prct_pev = e_rat;


    if verbose
       close all
       arrayfun(@(n) plot_rrrSummary(score{n},data(cur_fit).cvl_rrr{n},data(cur_fit).cvl_fa{n},bad_idx(n)),1:numel(bad_idx))

       %plot the distribution of model strength
       e_rat(bad_idx)=[]; withinSEM(bad_idx)=[]; 
       figure; hold on; 
       histogram(e_rat(withinSEM==1),'binwidth',0.01,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1]);
       histogram(e_rat(withinSEM==0),'binwidth',0.01,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8]);
       legend('Within SEM','Outside SEM','location','northeastoutside'); xlabel('performance fraction difference'); ylabel('# of subspaces')

       %confirm that strength is not just really weakly predictive ones
       temp = score; 
       temp(bad_idx)=[];
       temp = cellfun(@nanmean,temp); 
       figure; hold on;   
       histogram(temp(withinSEM==1),'binwidth',0.05,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.8 0.1 0.1]);
       histogram(temp(withinSEM==0),'binwidth',0.05,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor',[0.1 0.1 0.8]);
       legend('Within SEM','Outside SEM','location','northeastoutside'); xlabel('performance'); ylabel('# of subspaces')   
       
       %ratio versus predictive magnitude       
       figure; hold on; 
       lm = fitlm(mean_score(sum(bad_idx,2)==0),e_rat(sum(bad_idx,2)==0));
       plot(lm); title('fraction of full captured versus magnitude');
       
       DominantVersusPredictive(data(cur_fit).cvl_rrr, data(cur_fit).cvl_fa,sum(bad_idx_all{cur_fit}(:,[1:3,5:7]),2)>0,1);
    end
end %fit loop 



end %function end

%     for i = 1:numel(y) 
%         figure; hold on;
%         temp = abs(y{i});
%         histogram(temp); 
%         yvals = get(gca,'ylim');
%         plot([(prctile(temp,75)*5)+iqr(temp),(prctile(temp,75)*5)+iqr(temp)],yvals,'linestyle','--','color','r');
%         plot([(prctile(temp,25)*5)-iqr(temp),(prctile(temp,25)*5)-iqr(temp)],yvals,'linestyle','--','color','r');        
%         pause(); 
%         close
%     end




