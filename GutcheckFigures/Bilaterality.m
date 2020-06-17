function Bilaterality(motifs)

temp = [];
for i = 1:14
motif = reshape(squeeze(motifs(:,i,:)),[68, 68, size(motifs,3)]);
mL =(motif(:,1:34,:));
mR = (motif(:,35:end,:));
for j = 1:size(mR,3)
    mR(:,:,j) = fliplr(mR(:,:,j));
end

    


bad = sum([mL(:),mR(:)],2);

mL = mL(:);
mR = mR(:);
temp(i) = fisherZ(corr(mL(bad>0),mR(bad>0)));

end
rho_avg = fisherInverse(nanmean(temp));
rho_sem = fisherInverse(std(temp)./sqrt(numel(temp)));


