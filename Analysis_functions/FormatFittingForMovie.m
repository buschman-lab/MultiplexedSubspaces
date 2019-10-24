function FormatFittingForMovie(X,W,H,opts)
    %Smooth the original and the Ws
    X = conditionDffMat(X',opts.nanpxs{opts.block});
    X(isnan(X))=0;
    X = imgaussfilt3(X,[1, 1, 0.5]);
    X = reshape(X,size(X,1)^2,size(X,3));
    %Smooth and linearize W_all
    [nP, nK, nL] = size(W); %nK = #motifs nL = length of each motif
    nX = sqrt(nP);
    W_smooth = zeros(size(W));
    for cur_k = 1:nK
        temp = squeeze(W(:,cur_k,:));
        temp = reshape(temp,[nX nX nL]);
        temp = imgaussfilt3(temp,[1, 1, 0.5]);   
        temp = reshape(temp,[nP,nL]);
        W_smooth(:,cur_k,:) = temp;
    end  
    MakeExampleRepitoireReconstructionVideo(X,W_smooth,H,[],[opts.save_dir save_file_name '.avi'])
end