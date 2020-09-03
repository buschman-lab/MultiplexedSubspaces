function [total_entropy,windowed_entropy] = MotifEntropy(H,thresh,dur)
%Windowed Motif Entrop

[M,~] = size(H);

%threshold the H values with moving standard deviation
H_thresh=H;
for i =1:size(H,1)   
   temp = H(i,:);   
   if isempty(thresh) %use 2* std (e.g. if using Hs)
       threshval = 1*std(temp,[],2);
   else
       threshval = thresh; %use set value (e.g. if correlation). 
   end
   temp(temp<threshval)=0;
   H_thresh(i,:)=temp;  
end

% for computing the windowed entropy
[~,T] = size(H);

assert(floor(T/dur)==T/dur,'data not mutliple of window size');

num_chunks = floor(T/dur);


%get the maximum motif at each timepoint 
[~,X] = max(H_thresh,[],1);
bad_X = sum(H_thresh,1); %remove indices where no motifs are activity (otherwise that will be identified as motif 1 activity)
X(bad_X<=eps)=numel(unique(X))+1;

%reshape into chunks
X_chunked = reshape(X,[size(X,1),num_chunks,dur]);
% bad_X_chunked = reshape(bad_X,[size(X,1),num_chunks,dur]);

windowed_entropy = zeros(1,num_chunks);
for chunk = 1:size(X_chunked,2)
   Xtemp = squeeze(X_chunked(:,chunk,:));
%    Xtemp(squeeze(bad_X_chunked(:,chunk,:))<=eps) = [];
   X_uni = unique(Xtemp);
   X_uni_size = numel(X_uni);

   P = zeros(X_uni_size,1);

   for i = 1:X_uni_size
       P(i) = sum(Xtemp == X_uni(i));
   end

   P = P ./ numel(Xtemp);

   % Compute the Shannon's Entropy
   windowed_entropy(chunk) = -sum(P .* log2(P)); % 1.5
end


% total entropy
% bad_X = sum(H_thresh,1); %remove indices where no motifs are activity (otherwise that will be identified as motif 1 activity)
% X(bad_X==0)=[];

X_uni = unique(X);
X_uni_size = numel(X_uni);

P = zeros(X_uni_size,1);

for i = 1:X_uni_size
    P(i) = sum(X == X_uni(i));
end

P = P ./ numel(X);

% Compute the Shannon's Entropy
total_entropy = -sum(P .* log2(P)); % 1.5
 

    




end


