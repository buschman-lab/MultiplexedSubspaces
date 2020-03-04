function [X_recon, err, Wt, A, Ws] = NMFSpaceTime(X,varargin)
% X needs to be time by space
% Credit: 
% Camden MacDowell - timeless
%
% This function is based on the following work: 
% Onken et al., 2016 Using Matrix and Tensor Factorizations for the Single-Trial Analysis of
% Population Spike Trains
% and the seqNMF toolbox by Mackevicius and Bahle. 

% set adjustable parameters
params.maxiter = 250; 
params.verbose = 1; 
params.d_time = 1;
params.d_space = 1; 
params.err_tol = 0; 

params = ParseOptionalInputs(params,varargin);

N = size(X,2);
T = size(X,1);

%check non-negativity
assert(~any(X(:)<0),'Error: input data contains negative elements');

%initialize and normalize to frobenius norm 1 for either rows (Wt) or columns (Ws) 
Ws = normalize_frobenius(max(X(:))*rand(params.d_space,N),2); %row are pixels
Wt = normalize_frobenius(max(X(:))*rand(T,params.d_time),1); %columns are timepoints

%initialize the weights relating Ws and Wt
A=rand(params.d_time,params.d_space);  

%main loop 
err = NaN(1,params.maxiter);

%optionally run on the gpu
if isa(X,'gpuArray')
    Ws = gpuArray(Ws); 
    Wt = gpuArray(Wt); 
    A = gpuArray(A);  
    err = gpuArray(err);
end

for iter = 1:params.maxiter
    if mod(iter,10)==0 && params.verbose
        fprintf('\n\t fitting iteration %d',iter);        
    end
    %Check stopping criteria: e.g.
    
    %multiplicative update Ws
    G = Wt*A;
    Ws = Ws.*(G'*X)./((G'*G)*Ws+eps(Ws)); %add small numbers since can't update zeros
    
    %multiplicative update Wt
    V = A*Ws; 
    Wt = Wt.*(X*V')./(Wt*(V*V')+eps(Wt));        
    
    %update A
    A = A.*(Wt'*X*Ws')./((Wt'*Wt)*A*(Ws*Ws')+eps(A));
    
    %calculate error
    err(iter) = norm(X-Wt*A*Ws,'fro')^2;
    
    %normalize Ws and Wt
    [Ws, Ws_norm_mat] = normalize_frobenius(Ws,2);
    [Wt, Wt_norm_mat] = normalize_frobenius(Wt,1);
    
    %normalize A to maintain the error 
    for i = 1:size(Wt,2)
        for j = 1:size(Ws,1)
           A(i,j) = A(i,j).*(Wt_norm_mat(i)*Ws_norm_mat(j));
        end
    end
    
    %check end criteria
    if iter>1 && abs(err(iter)-err(iter-1))<params.err_tol
        break
    end

end %iter loop

if isa(X,'gpuArray')
    Ws = gather(Ws); 
    Wt = gather(Wt); 
    A = gather(A);  
    err = gather(err);
end

% remove not reached iterations
err(isnan(err)) =[];
X_recon = (Wt*A*Ws);

end %function
    



