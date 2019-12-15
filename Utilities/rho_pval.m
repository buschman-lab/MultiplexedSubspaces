function p = rho_pval(rho,N,tail)
if nargin <3; tail ='both'; end
%calculates pvalue from pearson's correlation coefficient using t
%distirbution

df = N-2;
tval = rho*sqrt(df)/sqrt(1-rho^2);

switch tail
    case 'right'
        p = tcdf(-tval, df);
    case 'left'
        p = tcdf(tval, df);
    case 'both'
        p = 2 * tcdf(-abs(tval), df);
end
       
end

