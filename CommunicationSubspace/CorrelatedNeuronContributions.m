function [rho,overlap,idx] = CorrelatedNeuronContributions(idx_dom,performance)
%Camden 
if nargin <2; performance = []; end

if ~isempty(performance) %order by performance
    [~,idx] = sort(nanmean(performance),'ascend');
    idx_dom = idx_dom(:,idx,:);
end
rng('default');
[n,m,~] = size(idx_dom);
ndim = 10; 
rho = NaN(m,m,n);
overlap = NaN(m,m,n);
for cur_rec = 1:n
    for cur_m  = 1:m
        for comp_m = 1:m
            if cur_m~=comp_m
                x = squeeze(idx_dom(cur_rec,cur_m,:));
                y = squeeze(idx_dom(cur_rec,comp_m,:));
                rho(cur_m,comp_m,cur_rec) = fisherZ(corr(x,y,'rows','complete','type','spearman'));
                overlap(cur_m,comp_m,cur_rec) = sum(ismember(x(1:ndim),y(1:ndim)))/ndim;   
            end
        end
    end
end

%to test the null here, you need to repeat the entire procedure 999 times,
%while you shuffle within each brain region. Then take the max correlation
%across all motifs for each shuffle as you permuted distribution

end %function 
