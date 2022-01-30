function [a_perm,b_perm] = neu_permutedCVs(x,y,num_cvs,n_perm)
%Camden MacDowell - timeless
%generally follows 'significantCV.m' ...
%returns 'num_cvs' CVs on data where the rows of x (i.e., neurons) are
%permuted. Used to determine significance of calculations computed on
%cv weightings

if nargin <4; n_perm = 1000; end

%create trail permuted of x
x_perm = NaN(size(x,2)*size(x,3),size(x,1),n_perm);
rng('default');
for i = 1:n_perm
   x_perm(:,:,i) = reshape(x(randperm(size(x,1),size(x,1)),:,:),[size(x,1),size(x,2)*size(x,3)])';
end

%concatentate across trials
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

%mean substracted
y = y-nanmean(y);
x_perm = x_perm-(nanmean(x_perm,1));

[a,b] = arrayfun(@(n) canoncorr(squeeze(x_perm(:,:,n)),y), 1:n_perm,'UniformOutput',0);


a = cellfun(@(x) x(:,1:num_cvs),a,'UniformOutput',0);
b = cellfun(@(x) x(:,1:num_cvs),b,'UniformOutput',0);
if num_cvs==1 %if only returning one CV then array
    a = cat(2,a{:});
    b = cat(2,b{:});
end

a_perm = a; 
b_perm = b; 

end %function end















