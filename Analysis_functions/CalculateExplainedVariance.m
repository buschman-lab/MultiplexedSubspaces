function [ExpVar_frame, ExpVar_all] = CalculateExplainedVariance(X,Xhat,Residuals)

%Calculate the variance
ExpVar_frame = zeros(1,size(Xhat,2));
for i = 1:size(Xhat,2)    
    ExpVar_frame(i) = 100*(1 - nanvar(Residuals(:,i))/nanvar(X(:,i)));
end

%Across all frames
ExpVar_all = 1 - nanvar(Residuals(:))./nanvar(X(:));

end