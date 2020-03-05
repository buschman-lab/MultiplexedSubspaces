function [stats_tab, anova_tbl] = StatByAnimal(x,id)
%Camden MacDowell - timeless
%Get the mean and CI of each statistic across animals. 

unique_id = unique(id);

%get the values for each animal
temp = arrayfun(@(y) x(id==y), unique_id, 'UniformOutput',0); 

clear stats
stats(1,:) = unique_id;
stats(2,:) = cellfun(@nanmedian,temp); %median
temp_ci = cellfun(@(y) bootci(1000,@nanmedian,y),temp,'UniformOutput',0); %median ci
stats(3,:) = cellfun(@(y) y(1),temp_ci); %median ci lower
stats(4,:) = cellfun(@(y) y(2),temp_ci); %median ci upper
stats(5,:) = cellfun(@nanmean,temp); %mean
temp_ci = cellfun(@(y) bootci(1000,@nanmean,y),temp,'UniformOutput',0); %mean ci
stats(6,:) = cellfun(@(y) y(1),temp_ci); %median ci lower
stats(7,:) = cellfun(@(y) y(2),temp_ci); %median ci upper

stats_tab = array2table(stats,'RowNames',{'ID','median','med_ci_lower','med_ci_higher','mean','mean_ci_lower','mean_ci_higher'});

%% run anova
[~,anova_tbl] = anova1(x,id,'off');


end %function end


