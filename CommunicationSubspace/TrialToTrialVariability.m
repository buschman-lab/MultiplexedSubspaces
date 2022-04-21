function ttt_var = TrialToTrialVariability(area_val,area_label,type,normtype)
%Camden - timeless
%get the trial to trial variabiltiy
if nargin <3; type = 'sum'; end
if nargin <4; normtype = 'mean'; end

ttt_var = NaN(1,numel(area_label));
for i = 1:numel(area_label)
    y = area_val{strcmp(area_label,area_label{i}),:};

    %normalize to baseline
    y = normalizeToBaseline(y,[1:2],normtype);

    %use post stimulus
    y = y(:,3:end,:);

    %subtract the psth
    y = y-nanmean(y,3);
    
    %flatten
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)]);    
        
    %variability across neurons    
    y = nanvar(y,[],2);
    
    switch type
        case 'mean'
            ttt_var(i) = nanmean(y);
        case 'sum'
            ttt_var(i) = nansum(y);
    end
    
end %subspace identification loop

end %function 


















