function loadings = computeLoadingPercentPower(V,W,H)

    K = size(H,1); 
    varv = sum(V(:).^2);
    loadings = zeros(1,K); 
    %Check For GPU Array
    if isa(V,'gpuArray'); loadings = gpuArray(loadings); end
    
    for fi = 1:K
        WH = helper.reconstruct(W(:,fi,:),H(fi,:)); 
        loadings(fi) = sum(2*V(:).*WH(:) - WH(:).^2)/varv;        
    end
    loadings(loadings<0)=0;
end