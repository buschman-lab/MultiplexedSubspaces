function LOOCVsubspace(x,y,verbose)
%Camden - timeless
%gets the strength of the loocv prediction of subspace y from subspace x

if nargin <3; verbose = 0; end

yhat = NaN(size(x,1),1);
for i = 1:size(x,1)
    idx = 1:size(x,1);
    train = idx(~ismember(idx,i));
    test = idx(ismember(idx,i));
    lm = fitlm(x(train),y(train));
    yhat(i) = predict(lm,x(test));
end

if verbose
   figure;
   lm = fitlm(y,yhat);
   plot(lm);
   ylabel('Predicted Weight');
   xlabel('True Weight');   
   title(sprintf('rsq%0.3f',lm.Rsquared.Ordinary))
end

end %function end
