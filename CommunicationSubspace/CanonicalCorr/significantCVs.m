function [a_full,b_full,U_full,V_full,rsq_xv,pval,rsq_xv_norm,aTheta,bTheta] = significantCVs(x,y,alpha,verbose)
%Camden MacDowell - timeless
%follows methods described in: https://www.nature.com/articles/s41467-020-17902-1#Sec13
%"single-trial cross-area neural population dynamics during long-term skill
%learning. Veuthey et al., Nat Comm. 2020. 
%performs cca and determine generalizability with cross validation and
%significant cannical variables with permutation
%inputs: 
%x and y are two tensors neurons x timepoits x trials. 
%alpha is the significance lvl

if nargin <3; alpha = 0.05; end 
if nargin <4; verbose = 1; end

%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%generate null distribution of y | trial-permuted
n_perm = 1000;
% y_perm = NaN(size(y,2)*size(y,3),10,n_perm);
y_perm = NaN(size(y,2)*size(y,3),size(y,1),n_perm);
rng('default');
for i = 1:n_perm
   %since pca performed over concatenated trials, this doesn't change the neuron coefficients, just the trial weightings
   y_perm(:,:,i) = reshape(y(:,:,randperm(size(y,3),size(y,3))),[size(y,1),size(y,2)*size(y,3)])';
end

%concatentate across trials and pca
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
% temp = U_full(:,1).*V_full(:,1);
% trueU = reshape(temp,[9,size(temp,1)/9]);

%run full cca
[a_full,b_full,r,U_full,V_full] = canoncorr(x,y); 


%get cross-validated R2 to test the generalization 
[rsq_xv,aTheta,bTheta] = xvalidateCCA(x,y,10);

%CCA on trial-permuted data for significance testing
% Option 1: conservative. no xvalidation for rsq_perm
[~,~,rsq_perm] = arrayfun(@(n) canoncorr(x,squeeze(y_perm(:,:,n))), 1:n_perm,'UniformOutput',0);
rsq_perm = cat(1,rsq_perm{:}).^2;
% Option 2: less conservative and slower. yes xvalidation for rsq_perm
% rsq_perm = arrayfun(@(n) xvalidateCCA(x,squeeze(y_perm(:,:,n)),10), 1:n_perm,'UniformOutput',0);
% rsq_perm = cat(1,rsq_perm{:});
%control for multiple comparisons by comparing all to the top perm CV
sig_thresh = prctile(rsq_perm(:,1),(1-alpha)*100);

if verbose %show the significant CVs
    figure('position',[681 764 423 215]); hold on; 
    histogram(rsq_perm(:,1),'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7)
    yyvals = get(gca,'ylim');
    arrayfun(@(xx) plot([xx xx], yyvals,'linewidth',0.5,'color',[0,0.2 0.82,0.4]),rsq_xv)
    plot([sig_thresh,sig_thresh],yyvals,'linewidth',2,'color','k','linestyle','--');
    ylabel('count');
    xlabel('Canonical Variable R^2','Interpreter','tex')
    title('Example Dataset CVs Significance','fontweight','normal');
end

%return significant cvs
sig_idx = rsq_xv>sig_thresh;
a_full = a_full(:,sig_idx);
b_full = b_full(:,sig_idx);
U_full = U_full(:,sig_idx);
V_full = V_full(:,sig_idx);
rsq_xv = rsq_xv(sig_idx);
rsq_xv_norm = rsq_xv-nanmedian(rsq_perm(:,1)); %how much larger relative to the trial permuted strength

%return a pvalue
pval = arrayfun(@(x) sum([rsq_perm(:,1)',x]>=x)/numel([rsq_perm(:,1)',x]), rsq_xv,'UniformOutput',1);

%output a warning if the first cv is less robust than subsequent
if ~isempty(rsq_xv) && sum(rsq_xv(1)<rsq_xv)>=1
    warning('first cv is less robust than subsequent');
end


end %function end















