function r = fisherInverse(z)


%Camden MacDowell 2022
f = @(x) ((exp(2*x))-1) / ((exp(2*x))+1);
r = arrayfun(@(x) (f(x)), z,'UniformOutput',1);

% %Camden MacDowell 2018
% r = zeros(size(z,1),size(z,2));
% for i = 1:size(z,1)
%     for j = 1:size(z,2)
%         r(i,j) = ((exp(2*z(i,j)))-1) / ((exp(2*z(i,j)))+1);
%     end
% end

end