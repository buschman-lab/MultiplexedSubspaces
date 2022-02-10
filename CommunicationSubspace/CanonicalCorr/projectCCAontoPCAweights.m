function a_full = projectCCAontoPCAweights(a,mu,coef)
%camden macdowell - timeless

a_full = a'*coef'+repmat(mu,size(a,2),1);

end %function end