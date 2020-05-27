function WindowedEntropy(H,dur)
%Windowed Motif Entrop

if nargin <2; dur =60*13; end 

dur =60*13; 
M = size(H,1); %number of motifs

windowed_data=WindowData(H,dur);


%entorpy across entire data
e_total = NaN(M,1);
for m = 1:M
   e_total(m) = entropy(H(m,:));
end

e = NaN(M,numel(windowed_data));
for i = 1:numel(windowed_data)
   for m = 1:M
      e(m,i) = entropy(windowed_data{i}(m,:));
   end  
end

%relative to entire recording
e_rel = e;
for m = 1:M
    e_rel(m,:) = smoothdata(e(m,:)/e_total(m),'movmean',4);
end
figure; hold on; 
plot(e_rel');
xlabel('time (min)');
ylabel('Entropy (relative to entire recording)')


%relative to entire recording
e_norm = e;
for m = 1:M
    e_norm(m,:) = smoothdata(e(m,:)-nanmean(e(m,:)),'movmean',3);
end

%figure 
plot(e_norm')

plot(e')



