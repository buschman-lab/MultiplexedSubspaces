function [X_recon,err,Wt,A,Ws, fit_to_data,opts] = DecomposeMotifs(varargin)
%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%for reproducibility
rng('default'); 

%initial parameters
fit_iter = 500; 
 
fprintf('\n\t Working on epoch %d \n',opts.block)

%Load the preprocessed data. Must be 1 x #rec cell arrays of data & nanpxs 
load([opts.bucket opts.base opts.data_file_name],'data');
X = data{opts.block};

%load the motifs
w = load([opts.bucket opts.base '/AnalyzedData_MesomappingManuscript_5_2019/TrainRepitoires/TrainingFit_Lambda4e-4_Kval28/TrainRepitoire_block_' num2str(opts.block)],'w');
w = w.w; 
d_motifs = size(w,2);

%% Fit to the raw dat using the same dimensions at cNMF
[fit_to_data.X_recon,  fit_to_data.err, fit_to_data.Wt, fit_to_data.A, fit_to_data.Ws] = NMFSpaceTime(X','maxiter',fit_iter,'verbose',0,'d_time',d_motifs,'d_space',d_motifs,'err_tol',0);
fit_to_data.pev = CalculateExplainedVariance(X,X-fit_to_data.X_recon');  

%% Fit to discovered motifs
%allign motifs to their maximal cross correlation with each other
w_alligned = AlignMotifs(w); 
X = reshape(w_alligned,size(w_alligned,1),size(w_alligned,2)*size(w_alligned,3));

% Main loop 
max_iter = 250; %iterations allowed to reach the explained variance of cnmf
d_time = 1; 
d_space = 1; 
pev_targ = 0.9;

%initital fit
dimensions = NaN(1,max_iter);
pev = NaN(1,max_iter);

for iter = 1:max_iter 
   if iter == 1      
       [A{1},Wt{1},Ws{1},~,err{1}] = sbtnmf(X,d_time,d_space,d_motifs,[],[],fit_iter,0,0,1);
       X_recon{1} = Reconstruct(A{1},Wt{1},Ws{1},d_motifs);
       dimensions(iter) = sum(A{1}(:)>eps);
       pev(iter) = CalculateExplainedVariance(X,X-X_recon{1}');
   else %greedily add dimensions
       %get change in pev if adding a space dimension
       [A{1},Wt{1},Ws{1},~,err{1}] = sbtnmf(X,d_time,d_space+1,d_motifs,[],[],fit_iter,0,0,1);
       X_recon{1} = Reconstruct(A{1},Wt{1},Ws{1},d_motifs);
       pev_temp{1} = CalculateExplainedVariance(X,X-X_recon{1}');

       %get change in pev if adding a time dimension
       [A{2},Wt{2},Ws{2},~,err{2}] = sbtnmf(X,d_time+1,d_space,d_motifs,[],[],fit_iter,0,0,1);
       X_recon{2} = Reconstruct(A{2},Wt{2},Ws{2},d_motifs);
       pev_temp{2} = CalculateExplainedVariance(X,X-X_recon{2}');
       
       %update 
       if pev_temp{1} >= pev_temp{2} %add space dimension       
           d_space = d_space+1;
           idx = 1; 
           fprintf('\n\t Iteration %d added space dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);                     
       else %add time dimension
           d_time = d_time+1;
           idx = 2;                 
           fprintf('\n\t Iteration %d added time dim. space = %d, time = %d PEV=%.4g',iter,d_space,d_time, pev_temp{1}*100);      
       end
       dimensions(iter) = sum(A{idx}(:)>eps); 
       pev(iter) = pev_temp{idx};
       
       if pev(iter)>=pev_targ %collected variables and end
           X_recon = X_recon{idx};
           err = err{idx};
           Wt = Wt{idx};
           A = A{idx};
           Ws =  Ws{idx};
           break
       elseif pev(iter)<pev_targ && iter == max_iter %populate with Nans if never reached the threshold
           X_recon = NaN;
           err = NaN;
           Wt = NaN;
           A = NaN;
           Ws =  NaN;
       else %keep iterating
       end
   end %if first iteration
end %iteration loop



end %end function loop

function X_recon = Reconstruct(A,Wt,Ws,d_motifs)
for i = 1:d_motifs
    temp = Wt*A(:,:,i)*Ws;
    if i ==1
        X_recon = temp;
    else
        X_recon = cat(1,X_recon,temp);
    end
end
end%function




