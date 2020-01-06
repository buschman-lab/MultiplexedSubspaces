function y = sem(x,dim)

if nargin <2; dim = 1; end

x(isnan(x))=[];
y = nanstd(x,[],dim)/sqrt(size(x,dim));
end

