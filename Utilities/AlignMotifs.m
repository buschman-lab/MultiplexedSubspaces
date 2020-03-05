function w_alligned = AlignMotifs(w)
%Camden MacDowell - timeless
% w is a P x K x L tensor of motifs
% pads and alligns motifs according to their cross correlation 

%for reproducibility
rng('default');

%pad motifs 
w_pad = cat(3, zeros(size(w)),w,zeros(size(w)));

%chose random template
template = squeeze(w_pad(:,randperm(size(w_pad,2),1),:));

%allign to best cross correlation with that template
w_alligned = w_pad; % allocate
for cur_k = 1:size(w,2)
    temp = squeeze(w_pad(:,cur_k,:));
    %find the best lag, including zero lag
    rho = NaN(1,size(w,3)+1);
    for cur_shift = 1:size(w,3)+1
        temp_shift = circshift(temp,cur_shift-1,2);
        rho(cur_shift) = corr(temp_shift(:),template(:));
    end
    [~, lag] = max(rho);
    %shift by the best lag
    w_alligned(:,cur_k,:) = circshift(temp,lag-1,2); %-1 since zero indexed    
end

%trim any timepoints with zero variance. 
temp = reshape(w_alligned,size(w_alligned,1)*size(w_alligned,2),size(w_alligned,3));
bad_tp = nanvar(temp,[],1)<=eps;

w_alligned(:,:,bad_tp')=[];

end %function end
