function [X_recon,err,Wt,A,Ws] = cNMFvsNMFSpaceTime(varargin)
%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%for reproducibility
rng('default'); 

fprintf('\n\t Working on epoch %d \n',opts.block)

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data');
X = data{opts.block};

%load the explained variance of the cNMF fit
w = load([opts.bucket opts.base '/AnalyzedData_MesomappingManuscript_5_2019/TrainRepitoires/TrainingFit_Lambda4e-4_Kval28/TrainRepitoire_block_' num2str(opts.block)],'w');

%% Restrict the dimensions to match cNMF
%initial parameters
max_iter = 100; %iterations allowed to reach the explained variance of cnmf
fit_iter = 500; 
d_time = size(w.w,3); 
d_space = size(w.w,3); 
[X_recon, ~,~, A, ~] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',d_time,'d_space',d_space,'err_tol',0);
pev = CalculateExplainedVariance(X,X-X_recon');  

%%


end %end function loop



%%

% %initital fit
% dimensions = NaN(1,max_iter);
% pev = NaN(1,max_iter);
% 
% for iter = 1:max_iter 
%    if iter == 1      
%         [X_recon{1}, ~,~, A{1}, ~] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',d_time,'d_space',d_space,'err_tol',0);
%         dimensions(iter) = sum(A{1}(:)>eps);
%         pev(iter) = CalculateExplainedVariance(X,X-X_recon{1}');        
%    else %greedily add dimensions
%        %get change in pev if adding a space dimension
%        [X_recon{1}, err{1}, Wt{1}, A{1}, Ws{1}] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',d_time,'d_space',d_space+1,'err_tol',0);
%        pev_temp{1} = CalculateExplainedVariance(X,X-X_recon{1}');
% 
%        %get change in pev if adding a time dimension
%        [X_recon{2}, err{2}, Wt{2}, A{2}, Ws{2}] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',d_time+1,'d_space',d_space,'err_tol',0);
%        pev_temp{2} = CalculateExplainedVariance(X,X-X_recon{2}');
%        
%        %update 
%        if pev_temp{1} >= pev_temp{2} %add space dimension       
%            d_space = d_space+1;
%            idx = 1; 
%            fprintf('\n\t Iteration %d added space dim. space = %d, time = %d PEV=%.4g',iter,d_time, d_space, pev_temp{1}*100);                     
%        else %add time dimension
%            d_time = d_time+1;
%            idx = 2;                 
%            fprintf('\n\t Iteration %d added time dim. space = %d, time = %d PEV=%.4g',iter,d_time, d_space, pev_temp{1}*100);      
%        end
%        dimensions(iter) = sum(A{idx}(:)>eps); 
%        pev(iter) = pev_temp{idx};
%        
%        if sum(d_time,d_space) == target_dim %collected variables and end
%            X_recon = X_recon{idx};
%            err = err{idx};
%            Wt = Wt{idx};
%            A = A{idx};
%            Ws =  Ws{idx};
%            break
%        end      
%    end %if first iteration
% end %iteration loop
% 
% %save off the data




