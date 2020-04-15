function X_recon = stNMF_Reconstruct(A,Wt,Ws)
for i = 1:size(A,3)
    temp = Wt*A(:,:,i)*Ws;
    if i ==1
        X_recon = temp;
    else
        X_recon = cat(1,X_recon,temp);
    end
end
end%function