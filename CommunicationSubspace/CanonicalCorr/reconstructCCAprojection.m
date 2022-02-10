function x = reconstructCCAprojection(a,u,mu,coef)
%camden macdowell - timeless

x_pca = u/a;
%To project from CCA to PCA to full space do the following
x = x_pca*coef'+repmat(mu,size(x_pca,1),1);

end %function end