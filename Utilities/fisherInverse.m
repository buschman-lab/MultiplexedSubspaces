function r = fisherInverse(z)

%Camden Macdowell 2018
%Note - z must be a vector
r = zeros(1,length(z));
for i = 1:length(z)
    r(i) = ((exp(2*z(i)))-1) / ((exp(2*z(i)))+1);
end

end