function CrossValidatedPCA(x,nfold,type)
%Camden - timeless
%x is the neuron x timepoint x trial matrix of trial to trial variability
%type is the type of cross validation timepoint selection
if nargin <3; type = 'trial'; end
if nargin <2; nfold = 10; end

switch type
    case 'trial'
        c = cvpartition(size(x,3),'KFold',nfold);
        for i = 1:nfold
            train = x(:,:,training(c,i));
            testdata = x(:,:,test(c,i));
            train = reshape(train,[size(train,1),size(train,2)*size(train,3)])';
            coef = pca(train);
            testdata = reshape(testdata,[size(testdata,1),size(testdata,2)*size(testdata,3)])';
            [~,scoretest,~,~,~,mu] = pca(testdata);
            %project test by training coefficients 
            yhat = arrayfun(@(n) 1-NormalizedSquaredError(testdata,scoretest(:,n)*coef(:,n)'+repmat(mu(:,n),size(scoretest,1),1)), 1:size(scoretest,2)); 
        end
    case 'random'
        c = cvpartition(size(x,3)*size(x,2),'KFold',nfold);
    otherwise 
        error('unknown type');
end

for i = 1:nfold
   training = 

end

training(c,2)-training(c,1)








end
