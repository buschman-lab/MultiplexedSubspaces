function [explained, avg_trial] = TrialPEV(x,trial_idx)
%camden macdowell 2019
%compute the percent explained variance of each trial type to the total
%variance of the x. 
%x is a (t x 1) vector. If x is a t x n matrix. then trial based explained
%variance is performed per column. 

unique_trials = unique(trial_idx);
temp_x = x;
avg_trial = NaN(numel(unique_trials),size(temp_x,2));
for trial_type = 1:numel(unique_trials)    
    avg_trial(trial_type,:) = nanmean(x(trial_idx==unique_trials(trial_type),:),1);
    temp_x(trial_idx==unique_trials(trial_type),:) = temp_x(trial_idx==unique_trials(trial_type),:)-avg_trial(trial_type,:);
end
explained = 1-nanvar(temp_x)./nanvar(x);

end
