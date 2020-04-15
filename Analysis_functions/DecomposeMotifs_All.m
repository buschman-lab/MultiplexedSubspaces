function DecomposeMotifs_All(varargin)
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
end

%Set options
opts = SetAnalysisOptions();
opts = ParseOptionalInputs(opts,varargin); 

%Save off
save_dir = [opts.bucket opts.base '/TrainRepitoires/SpaceTimeNMF_AllMotifs/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%for reproducibility
rng('default'); 

%initial parameters
fit_iter = 250; 


%Load the palligned motifs
w_alligned = load([opts.bucket opts.base '/AnalyzedData_MesomappingManuscript_5_2019/TrainRepitoires/TrainingFit_Lambda4e-4_Kval28/AverageDPs_suppliment_1.mat'],'W_smooth_alligned');
idx = load([opts.bucket opts.base '/AnalyzedData_MesomappingManuscript_5_2019/TrainRepitoires/TrainingFit_Lambda4e-4_Kval28/ClusteredDPs_Paramset.mat'],'idx_louvain');
idx = idx.idx_louvain(:,2);
w_alligned =w_alligned.W_smooth_alligned;

unique_m = unique(idx);
subset = [];
for i = 1:14
    temp = find(idx == unique_m(i));
    temp = temp(randperm(numel(temp),ceil(numel(temp)/4)));
    subset = cat(1,subset, temp);
end
w_alligned = w_alligned(:,subset,:);   
   
d_motifs = size(w_alligned,2);

%% Fit to discovered motifs
X = reshape(w_alligned,size(w_alligned,1),size(w_alligned,2)*size(w_alligned,3));

bad_dim = isnan(nanvar(X,[],2)) | nanvar(X,[],2)<=eps;
X(bad_dim,:) = [];
X(isnan(X))=0;

% Main loop 
max_iter = 75; %iterations allowed to reach the explained variance of cnmf
d_time = 1; 
d_space = 1; 
pev_targ = 0.95;

%initital fit
pev = NaN(1,max_iter);
for iter = 1:max_iter    
   if iter == 1    
       fprintf('\nWorkingOnIteration1');
       [A{1},Wt{1},Ws{1},~,err{1}] = sbtnmf(X,d_time,d_space,d_motifs,[],[],fit_iter,0,0,1);
       X_recon{1} = Reconstruct(A{1},Wt{1},Ws{1},d_motifs);
       pev(iter) = CalculateExplainedVariance(X,X-X_recon{1}');
       idx = 1;

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
       pev(iter) = pev_temp{idx};
   end
       
   %Save off
   X_recon_temp = X_recon{idx};
   err_temp = err{idx};
   Wt_temp = Wt{idx};
   A_temp = A{idx};
   Ws_temp =  Ws{idx};
   filename = [save_dir sprintf('iter_%d.mat',iter) '.mat'];
   save(filename,'X_recon_temp','err_temp','Wt_temp','A_temp','Ws_temp','pev','iter','-v7.3'); 
   if pev(iter)>=pev_targ %collected variables and end
       break
   end
   
end %iteration loop



end %end function loop
% 
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




