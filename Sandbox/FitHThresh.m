function FitHThresh(h,motifs,data)

if nargin <2
    thresh = (0:0.25:4);
end

pev = zeros(numel(thresh),size(motifs,2));
%sweep threshold level
for i = 1:numel(thresh)
    fprintf('\n\t Working on threshold level %d or %d',i,numel(thresh));
    for j = 1:size(motifs,2) 
        
        X = tensor_convolve(motifs(:,j,:),h(j,:));
        [pks,locs] = findpeaks(h(j,:),'Threshold',thresh(i)*nanstd(h(j,:)),'MinPeakDistance',1*15);
        temp = zeros(size(h(j,:)));
        temp(locs)=pks;
        Xhat = tensor_convolve(motifs(:,j,:),temp);
        pev(i,j) = CalculateExplainedVariance(X,X-Xhat);
    end
end


end

figure; hold on;
plot(h(j,:),'linewidth',2);
plot(locs,pks,'*');
