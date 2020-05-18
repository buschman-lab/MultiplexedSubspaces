function [w, h, loss] = fpCNMF(X,varargin)

%parse optional inputs
opts.K = 5;
opts.L = 20;
opts = ParseOptionalInputs(opts,varargin);

%initialize
[N,T] = size(X);
w_init = max(X(:))*rand(N, opts.K, opts.L);%max(X(:))*rand(N, K, L);
h_init = max(X(:))*rand(opts.K,T)./(sqrt(T/3)); %max(X(:))*rand

tic
L = opts.L;
K = opts.K;
%Run CNMF_ALS + CNMF_pMU
for iter = 1:20 %% add the padding
    if iter ==1
        X_pad = [zeros(N,L),X,zeros(N,L)];
        w = w_init; 
        h = [zeros(K,L),h_init,zeros(K,L)];
        Xhat = helper.reconstruct(w,h);       
        cost = sqrt(mean((X_pad(:)-Xhat(:)).^2));
        [w,h] = CNMF_ALS(X_pad,w,h);   
    else
        [w,h] = CNMF_ALS(X_pad,w,h);   
    end
       
    Xhat = helper.reconstruct(w,h);   
    cost(iter+1) = sqrt(mean((X_pad(:)-Xhat(:)).^2));
%     loss(iter+1) = norm(X_pad-Xhat,'fro')/norm(X,'fro');
    %fit until the loss plateuas    
end
h = h(:,L+1:end-L);
t = toc;

%switch the CNMF_pMU 
[w, h, temp] = seqNMF(X, ...    % X is the data matrix
  'K', opts.K, 'L', opts.L, 'lambda', 0, ...        % Other inputs optional
  'showPlot', 0, 'maxiter', 10, 'tolerance', 0, 'shift', 0, ... 
  'lambdaL1W', 0, 'lambdaL1H', 0, ...
  'lambdaOrthoH', 0, 'lambdaOrthoW', 0,...
  'W_init',w,'H_init',h);
temp = temp(2:end);


%Run CNMF_MU_Penalized

%switch the CNMF_pMU 
tic
[w, h, cost2] = seqNMF(X, ...    % X is the data matrix
  'K', opts.K, 'L', opts.L, 'lambda', 0, ...        % Other inputs optional
  'showPlot', 0, 'maxiter', 20, 'tolerance', 0, 'shift', 0, ... 
  'lambdaL1W', 0, 'lambdaL1H', 0, ...
  'W_init', w_init,'H_init', h_init,...
  'lambdaOrthoH', 0, 'lambdaOrthoW', 0);
cost2 = cost2(1:end);
t2 = toc;
  

figure; cost = [cost, temp'];
plot(cost); hold on;
plot(cost2); hold on;



figure;hold on; plot((0:20)*t1,cost,'k'); plot((0:20)*t2,cost2,'b');set(gca,'yscale','log'); ylim([0 3])
plot(cost); hold on;
plot(cost2); hold on;

