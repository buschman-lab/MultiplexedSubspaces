function y = sem(x,dim)

if nargin <2; dim = 1; end

y = nanstd(x,[],dim)./sqrt(sum(~isnan(x),dim));
end

