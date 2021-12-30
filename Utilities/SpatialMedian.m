function [X] = SpatialMedian(X, kernelsize)
%This MATLAB function performs median filtering of the matrix A in two dimensions
%while *ignoring* NaNs (based on discussions here
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/251787)

%A is a 3d matrix, sz is tile size, M is filtered A.

%novariance = NaN
if sum(isnan(X(:)))==0 && size(X,3)>1
    temp = reshape(X,[size(X,1)*size(X,2),size(X,3)]);
    temp(nanvar(temp,[],2)<=eps,:) = NaN;
    X = reshape(temp,[size(X,1),size(X,2),size(X,3)]);    
end

%get nan values; 
mask = isnan(X);

for i = 1:size(X,3)
    if nargin<2
        kernelsize = 5;
    end
    if length(kernelsize)==1
        kernelsize = [kernelsize kernelsize];
    end

    if any(mod(kernelsize,2)==0)
        error('kernel size SZ must be odd)')
    end
    A = X(:,:,i);
    margin=(kernelsize-1)/2;
    AA = nan(size(A)+2*margin);
    AA(1+margin(1):end-margin(1),1+margin(2):end-margin(2))=A;
    [iB,jB]=ndgrid(1:kernelsize(1),1:kernelsize(2));
    is=sub2ind(size(AA),iB,jB);
    [iA, jA]=ndgrid(1:size(A,1),1:size(A,2));
    iA=sub2ind(size(AA),iA,jA);
    idx=bsxfun(@plus,iA(:).',is(:)-1);

    B = sort(AA(idx),1);
    j=any(isnan(B),1);
    last = zeros(1,size(B,2))+size(B,1);
    [~, last(j)]=max(isnan(B(:,j)),[],1);
    last(j)=last(j)-1;

    M = nan(1,size(B,2));
    valid = find(last>0);
    mid = (1 + last)/2;
    i1 = floor(mid(valid));
    i2 = ceil(mid(valid));
    i1 = sub2ind(size(B),i1,valid);
    i2 = sub2ind(size(B),i2,valid);
    M(valid) = 0.5*(B(i1) + B(i2));
    X(:,:,i) = reshape(M,size(A));
end
%return normal masked
X(mask)=NaN;
end % medianna

