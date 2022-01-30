function GutCheckSimilarityInMotifOnsets(data)
%gut check that motif onsets are different
%data variable from ANalyze CCs. could also just use motifOnsets
onset = data(1).motif_onset;
n = max(cellfun(@(x) max(x), onset));
temp = zeros(numel(onset),n+16);
for i = 1:numel(onset)
    for j = 1:numel(onset{i})
        temp(i,onset{i}(j):(onset{i}(j)+15))=1;
    end
end
figure('units','normalized','position',[0.0063 0.5657 0.9865 0.3398]); 
imagesc(temp); colormap(flipud(gray))
set(gca,'xtick',1:(10*60*15):size(temp,2))
set(gca,'XTickLabel',round(get(gca,'xtick')/(60*15)))
xlabel('time (min)'); ylabel('motif')
title('Motif onsets throughout recording','FontWeight','normal')

%% zoom in 
figure('units','normalized','position',[0.0063 0.5657 0.9865 0.3398]); 
imagesc(temp); colormap(flipud(gray))
set(gca,'xtick',1:(1*60*15):size(temp,2))
set(gca,'XTickLabel',round(get(gca,'xtick')/(60*15)))
xlabel('time (min)'); ylabel('motif')
title('(Zoomed) Motif onsets throughout recording','FontWeight','normal')
%zoom in to minutes 8.5-9.5
xlim([8.5*60*15, 10.5*60*15]);

%% get the amount of overlap between motifs
overlap_mat = NaN(numel(onset),numel(onset));
for i = 1:numel(onset)
    for j = 1:numel(onset)        
       m1 = temp(i,:); 
       m2 = temp(j,:); 
       %get the percentage of m1 active that is shared with m2       
       overlap_mat(i,j) = 100*(sum(sum(cat(1,m1,m2),1)==2)/sum(m1));
    end
end
overlap_mat(1:1+size(overlap_mat,1):end) = NaN;
figure; hold on; 
imagesc(overlap_mat,[0 50]); colormap(flipud(gray))
colorbar; axis square;
xlim([0.5 numel(onset)+0.5]);
ylim([0.5 numel(onset)+0.5]);
xlabel('motif'); ylabel('motif')
title(sprintf('%% overlap of motifs (avg=%0.2f%%)',nanmean(overlap_mat(:))),'fontweight','normal')

end