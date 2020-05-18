%testingCode

number_of_seqences = 5;
T = 5000; % length of data to generate
Nneurons = 150*ones(number_of_seqences,1); % number of neurons in each sequence
Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
NeuronNoise = 0.001; % probability of added noise in each bin
SeqNoiseTime = zeros(number_of_seqences,1); % Jitter parameter = 0%
SeqNoiseNeuron = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
data = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,0,0,0,0,0);

%initialize W and H

[N,T] = size(data);
K = 3;
L = 11;

% w = rand(N,K,L);
% h = rand(K,T);
% est = helper.reconstruct(w,h);
% alpha = (reshape(data,N*T,1)'* reshape(est,N*T,1))/norm(est)^2;
% w = w*sqrt(abs(alpha));
% h = h*sqrt(abs(alpha));

[w, h] = seqNMF(data, ...    % X is the data matrix
      'K', 3, 'L', 11, 'lambda', 0, ...        % Other inputs optional
      'showPlot', 1, 'maxiter', 20, 'tolerance', 0, 'shift', 0, ... 
      'lambdaL1W', 0, 'lambdaL1H', 0, ...
      'lambdaOrthoH', 0, 'lambdaOrthoW', 0);


%%


[N,T] = size(data);
K = 7;
L = 4;
w = max(data(:))*rand(N, K, L);%max(X(:))*rand(N, K, L);
h = max(data(:))*rand(K,T)./(sqrt(T/3)); %max(X(:))*rand

%update W; 
[N,K,L] = size(w);

Xhat = helper.reconstruct(w,h);
cost = sqrt(mean((data(:)-Xhat(:)).^2));
resids = data - helper.reconstruct(w,h);
loss(1) = norm(resids,'fro')/norm(data,'fro');
%Now update H (one at a time)

for iter = 1:25
    %make H into a block matrix format
    h_block = BlockH(h,L);  

    %Compute the nnls
%     if iter ==1
        w_stacked = nnlsm_blockpivot(h_block', data',0);
%     else
%         w_stacked = nnlsm_blockpivot(h_block', data',0,w_stacked); 
%     end
    
    %Unblock w
    w = ConditionW(w_stacked',[N,K,L]); 

    for t = 1:T 
        last = min(t+L-1,T);
        block_size = last - t + 1;
               
        b = data(:,t:last);
        b = b(:);
        w_blocked = ConditionW(w,[N*L,K]);
%         if iter ==1
        h(:,t) = nnlsm_blockpivot(w_blocked(1:block_size*N,:), b,0);        
%         else
%             h(:,t) = nnlsm_blockpivot(w_blocked(1:block_size*N,:), b,0,h(:,t));        
%         end
    end

    loss(iter+1) = norm(data-helper.reconstruct(w,h),'fro')/norm(data,'fro');
end

%%


%     %Now update H (one at a time)
%     resids = data - helper.reconstruct(w,h);
%     for t = 1:T
%         last = min(t+L-1,T);
%         block_size = last - t + 1;
%         
%         %remove the contribution to the residual
%         for k = 1:K
%            resids(:,t:last) = resids(:,t:last) - h(k,t)*squeeze(w(:,k,1:block_size));
%         end    
%         w_block = ConditionW(w,[N*L,K]);
%         b = resids(:,t:last); b = b(:);
%         %update single column of h
%         h(:,t) = nnlsm_blockpivot(w_block(1:block_size*N,:), -b);
% 
%         %update residuals
%         for k = 1:K
%            resids(:,t:last) = resids(:,t:last) + h(k,t)*squeeze(w(:,k,1:block_size));
%         end  
% 
%     end

% 
%     for cur_L = 1:L
%         inds = cur_L:L:T-L+1;
% 
%         %remove the contribution to the residual
%         for k = 1:K
%             for t = inds
%                 resids(:,t:t+L-1) = resids(:,t:t+L-1)-h(k,t)*squeeze(w(:,k,:));
%             end           
%         end    
%         w_block = ConditionW(w,[N*L,K]);
%         
%         b = zeros(N*L,numel(inds));
%         for i = 1:numel(inds)
%             t = inds(i);
%             temp = resids(:,t:t+L-1);
%             b(:,i)=temp(:);
%         end
%                
%         %update block of h
%         h(:,inds) = nnlsm_blockpivot(w_block, -b);
% 
%         %update residuals
%         for k = 1:K
%             for t = inds
%                 resids(:,t:t+L-1) = resids(:,t:t+L-1)+h(k,t)*squeeze(w(:,k,:));
%             end           
%         end         
%     end












