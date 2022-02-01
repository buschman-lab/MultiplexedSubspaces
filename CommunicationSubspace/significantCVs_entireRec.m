function [a_full,b_full,U_full,V_full,r,r_norm] = significantCVs_entireRec(x,y,alpha)
%Camden MacDowell - timeless
%x and y is the n (neuron) x t (time) activity of two neural populations 

if nargin <3; alpha = 0.05; end 

%transpose a mean substracted
x = x'; y = y';
x = x-nanmean(x);
y = y-nanmean(y);

%run full cca
[a_full,b_full,r,U_full,V_full,stats] = canoncorr(x,y);    

%get rel strength of cv is on performed on the residual of the previous 
r = r.^2;
r_norm = r/sum(r);

rng('default')
[~,~,r_perm] = arrayfun(@(n) canoncorr(x(randperm(size(x,1),size(x,1)),:),y),1,'UniformOutput',0);    


%return significant cvs
sig_idx = stats.p<alpha;
a_full = a_full(:,sig_idx);
b_full = b_full(:,sig_idx);
U_full = U_full(:,sig_idx);
V_full = V_full(:,sig_idx);
r_norm = r_norm(sig_idx);
r = r(sig_idx);



end %function end















