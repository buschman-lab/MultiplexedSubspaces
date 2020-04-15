function [fit_to_data,fit_full,opts] = DecomposeMotifs(varargin)

% [X_recon,err,Wt,A,Ws, fit_to_data,opts,dimensions,pev,rank_a] = old
% version of the code
%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%for reproducibility
rng('default'); 

%initial parameters
fit_iter = 1500; 
 
fprintf('\n\t Working on epoch %d \n',opts.block)

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data');
X = data{opts.block};

%load the motifs
temp = load([opts.bucket opts.base '/AnalyzedData_MesomappingManuscript_5_2019/TrainRepitoires/TrainingFit_Lambda4e-4_Kval28/TrainRepitoire_block_' num2str(opts.block)],'w','ExpVar_all');
pev_targ = temp.ExpVar_all;
w = temp.w; 

d_motifs = size(w,2);

%% Fit to the raw dat using the same dimensions at cNMF
% [fit_to_data.X_recon,  fit_to_data.err, fit_to_data.Wt, fit_to_data.A, fit_to_data.Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',1,'d_time',d_motifs,'d_space',d_motifs,'err_tol',0);
% fit_to_data.pev = CalculateExplainedVariance(X,X-fit_to_data.X_recon');  
fit_to_data = [];

%% Iteratively fit
max_dim = 251;
fit_full = struct();
iter = 1;
for cur_dim = [1:19,20:10:max_dim] %subsample    
   fprintf('\nWorking on dimension %d of %d',cur_dim,max_dim);
   [fit_full(iter).X_recon,  fit_full(iter).err, fit_full(iter).Wt, fit_full(iter).A, fit_full(iter).Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',cur_dim,'d_space',cur_dim,'err_tol',0);
   fit_full(iter).pev = CalculateExplainedVariance(X,X-fit_full(iter).X_recon'); 
   fit_full(iter).dim = cur_dim; 
   iter = iter+1;
end

%%
% [fit_full.X_recon,  fit_full.err, fit_full.Wt, fit_full.A, fit_full.Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',1,'d_time',size(X,2),'d_space',size(X,1),'err_tol',0);
% fit_full.pev = CalculateExplainedVariance(X,X-fit_full.X_recon');  

% %% Greedy Fit
% max_iter = 100;
% d_time = 1; 
% d_space = 1; 
% pev = NaN(1,max_iter);
% dimensions = NaN(max_iter,2);
% for iter = 1:max_iter
%      if iter == 1   
%         [X_recon,  ~, Wt, A, Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',1,'d_time',d_time,'d_space',d_space,'err_tol',0);
%         pev(iter) = CalculateExplainedVariance(X,X-fit_full.X_recon');    
%      else %greedily add dimensions
%         %get change in pev if adding a space dimension
%         [X_recon,  ~, Wt, A, Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',1,'d_time',d_time,'d_space',d_space+1,'err_tol',0);
%         pev_temp{1} = CalculateExplainedVariance(X,X-X_recon{1}');
% 
%         [X_recon,  ~, Wt, A, Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',1,'d_time',d_time+1,'d_space',d_space,'err_tol',0);
%         pev_temp{2} = CalculateExplainedVariance(X,X-X_recon{1}');
%        
%         %update 
%         if pev_temp{1} >= pev_temp{2} %add space dimension       
%            d_space = d_space+1;
%            idx = 1; 
%            fprintf('\n\t Iteration %d added space dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);                     
%         else %add time dimension
%            d_time = d_time+1;
%            idx = 2;                 
%            fprintf('\n\t Iteration %d added time dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);      
%         end
%         pev(iter) = pev_temp{idx};
%      end
%      dimensions(iter,:) = [d_space,d_time];
%       
%      if pev(iter)>=pev_targ %collected variables and end
%          break
%      end
%        
%    end %if first iteration
% end %iteration loop

%% THIS IS OLDER CODE THAT DECOMPOSES THE MOTIFS INTO SPACE AND TIME AS WELL. 
% %% Fit to discovered motifs
% %allign motifs to their maximal cross correlation with each other
% w_alligned = AlignMotifs(w); 
% X = reshape(w_alligned,size(w_alligned,1),size(w_alligned,2)*size(w_alligned,3));
% 
% % Main loop 
% max_iter = 1; %iterations allowed to reach the explained variance of cnmf
% d_time = 1; 
% d_space = 1; 
% pev_targ = 0.9;
% 
% %initital fit
% dimensions = NaN(3,max_iter);
% pev = NaN(1,max_iter);
% rank_a = cell(1,max_iter);
% for iter = 1:max_iter 
%    if iter == 1      
%        [A{1},Wt{1},Ws{1},~,err{1}] = sbtnmf(X,d_time,d_space,d_motifs,[],[],fit_iter,0,0,1);
%        X_recon{1} = Reconstruct(A{1},Wt{1},Ws{1},d_motifs);
%        rank_a{iter}(1,:) = arrayfun(@(n) rank(squeeze(A{1}(:,:,n))),1:size(A{1},3),'UniformOutput',1);
%        rank_a{iter}(2,:) = arrayfun(@(n) size(squeeze(A{1}(:,:,n)),1),1:size(A{1},3),'UniformOutput',1);
%        rank_a{iter}(3,:) = arrayfun(@(n) size(squeeze(A{1}(:,:,n)),2),1:size(A{1},3),'UniformOutput',1);
%        dimensions(:,iter) = cat(1,sum(A{1}(:)>eps),d_time,d_space);
%        pev(iter) = CalculateExplainedVariance(X,X-X_recon{1}');
%        
%        
%    else %greedily add dimensions
%        %get change in pev if adding a space dimension
%        [A{1},Wt{1},Ws{1},~,err{1}] = sbtnmf(X,d_time,d_space+1,d_motifs,[],[],fit_iter,0,0,1);
%        X_recon{1} = Reconstruct(A{1},Wt{1},Ws{1},d_motifs);
%        pev_temp{1} = CalculateExplainedVariance(X,X-X_recon{1}');
% 
%        %get change in pev if adding a time dimension
%        [A{2},Wt{2},Ws{2},~,err{2}] = sbtnmf(X,d_time+1,d_space,d_motifs,[],[],fit_iter,0,0,1);
%        X_recon{2} = Reconstruct(A{2},Wt{2},Ws{2},d_motifs);
%        pev_temp{2} = CalculateExplainedVariance(X,X-X_recon{2}');
%        
%        %update 
%        if pev_temp{1} >= pev_temp{2} %add space dimension       
%            d_space = d_space+1;
%            idx = 1; 
%            fprintf('\n\t Iteration %d added space dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);                     
%        else %add time dimension
%            d_time = d_time+1;
%            idx = 2;                 
%            fprintf('\n\t Iteration %d added time dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);      
%        end
%        rank_a{iter}(1,:) = arrayfun(@(n) rank(squeeze(A{idx}(:,:,n))),1:size(A{idx},3),'UniformOutput',1);
%        rank_a{iter}(2,:) = arrayfun(@(n) size(squeeze(A{idx}(:,:,n)),1),1:size(A{idx},3),'UniformOutput',1);
%        rank_a{iter}(3,:) = arrayfun(@(n) size(squeeze(A{idx}(:,:,n)),2),1:size(A{idx},3),'UniformOutput',1);
%        dimensions(:,iter) = cat(1,sum(A{idx}(:)>eps),d_time,d_space);     
%        pev(iter) = pev_temp{idx};
%        
%        
%        if pev(iter)>=pev_targ %collected variables and end
%            X_recon = X_recon{idx};
%            err = err{idx};
%            Wt = Wt{idx};
%            A = A{idx};
%            Ws =  Ws{idx};
%            break
%        elseif pev(iter)<pev_targ && iter == max_iter %populate with Nans if never reached the threshold
%            X_recon = NaN;
%            err = NaN;
%            Wt = NaN;
%            A = NaN;
%            Ws =  NaN;
%        else %keep iterating
%        end
%    end %if first iteration
% end %iteration loop



end %end function loop
% 
% function X_recon = Reconstruct(A,Wt,Ws,d_motifs)
% for i = 1:d_motifs
%     temp = Wt*A(:,:,i)*Ws;
%     if i ==1
%         X_recon = temp;
%     else
%         X_recon = cat(1,X_recon,temp);
%     end
% end
% end%function




