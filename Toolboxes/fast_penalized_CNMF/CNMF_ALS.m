function [w,h] = CNMF_ALS(data,w,h)
%update w (N x K x L tensor) and h (K x T matrix) using ALS. Performs a
%single iteration.

[N, K, L] = size(w);
[~, T] = size(h);

%make H into a block matrix format
h_block = BlockH(h,L);  

%Compute the nnls
w_stacked = nnlsm_blockpivot(h_block', data',0); 

%stack w
w = ConditionW(w_stacked',[N,K,L]); 

for t = 1:T 
    last = min(t+L-1,T);
    block_size = last - t + 1;

    b = data(:,t:last);
    b = b(:);
    w_blocked = ConditionW(w,[N*L,K]);
    h(:,t) = nnlsm_blockpivot(w_blocked(1:block_size*N,:), b,0);        
end

% loss = norm(data-helper.reconstruct(w,h),'fro')/norm(data,'fro');

end

