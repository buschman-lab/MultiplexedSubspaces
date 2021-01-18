function [p, mean_perm, mean_true] = OneSamplePermuationTest(x,y,tail)

rng('default')
true_diff = x-y;
mean_true = nanmean(true_diff);
mean_perm = NaN(1,1000);
for i = 1:1000
    idx = datasample([-1,1],N)'; %randomly sample from {-1, 1}
    mean_perm(i) = nanmean(true_diff.*idx);
end

switch tail
    case 'right'
        p = sum(mean_perm>=mean_true)/(1000+1);
    case 'left'       
        p = sum(mean_perm<=mean_true)/(1000+1);
end


end
