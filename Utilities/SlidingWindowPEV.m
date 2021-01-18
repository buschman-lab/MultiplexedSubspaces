function [pev,var_level] = SlidingWindowPEV(X,Xhat,dur,shift)
    [~,T] = size(X);
    
%     chunks = cat(1,1:shift:T-shift+1,shift:shift:T);
    chunks = cat(1,1:shift:T-dur+1,(1:shift:T-dur+1)+dur-1);

    pev = NaN(1,size(chunks,2));
    var_level= NaN(1,size(chunks,2));
    for i = 1:size(chunks,2)
       temp_X = X(:,chunks(1,i):chunks(2,i));
       temp_Xhat = Xhat(:,chunks(1,i):chunks(2,i));
       pev(i) = CalculateExplainedVariance(temp_X,temp_X-temp_Xhat);
       var_level(i) = nanvar(temp_X(:));
    end
        
end