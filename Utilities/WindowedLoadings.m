function loadings = WindowedLoadings(w,H,data,dur)
    [~,T] = size(data);

    assert(floor(T/dur)==T/dur,'data not mutliple of window size');

    num_chunks = floor(T/dur);
    if mod(num_chunks,2)% need even number to equally split test/train
        num_chunks = num_chunks-1; 
    end

    %reshape into chunks
    X = reshape(data,[size(data,1),num_chunks,dur]);
    H_chunked = reshape(H,[size(H,1),num_chunks,dur]);

    %compute the loadings for each chunk
    K = size(H,1); 
    loadings = zeros(K,num_chunks);
    for chunk = 1:num_chunks
       %loadings of each w   
       temploadings = zeros(1,K);   
       for i = 1:K
           Xhat = tensor_convolve(w(:,i,:),squeeze(H_chunked(i,chunk,:))');
           temploadings(i) = CalculateExplainedVariance(squeeze(X(:,chunk,:)),squeeze(X(:,chunk,:))-Xhat);     
       end
       loadings(:,chunk) = temploadings/sum(temploadings);
    end
end