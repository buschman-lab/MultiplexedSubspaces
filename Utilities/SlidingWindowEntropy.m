function [ent_motifs, ent_transitions, num_transitions] = SlidingWindowEntropy(H,nstd,dur,shift)
    %entropy of motifs is the relative entropy (KL-divergence) computed on the probability
    %of each motif being 'dominant' for each point in time. For example, 
    %If there is an equal probabily of all motifs occuring then ent_motifs
    %will be high. If the probability is skewed then the entropy will be
    %lower. ent_transition is the relative entropy  entropy computed using the
    %probability of each transition between motifs. For example, if a small
    %group of transitions dominates the transition matrix, then the entropy
    %will be low. If it's more disributed probabiltiy transitions between 
    %motifs then entropy will be high. 

    [N,T] = size(H);
    
    %split into windows
    chunks = cat(1,1:shift:T-dur+1,(1:shift:T-dur+1)+dur-1);

    %threshold the H's
    H_thresh = H;
%     H_thresh(H<repmat(nstd*std(H,[],2),1,T))=0;
    H_thresh(H<(nstd*movstd(H,dur,0,2)))=0;
    
    %get the dominant motif at each timepoint
    [~,idx_dom] = max(H_thresh,[],1);
    idx_dom(sum(H_thresh,1)==0)=N+1; %Add a 'null state' for when no motifs are active
    
   
    ent_motifs = NaN(1,size(chunks,2));
    ent_transitions = NaN(1,size(chunks,2));
    num_transitions = NaN(1,size(chunks,2));
    for i = 1:size(chunks,2)
       temp_idx = idx_dom(:,chunks(1,i):chunks(2,i)); 
       
       %compute the entropy of motifs
       p_m = NaN(1,N+1);
       for j = 1:N+1    
           p_m(j) = (sum(temp_idx==j)+eps)/numel(temp_idx); 
       end       
       ent_motifs(i) = sum(p_m .* (log2(p_m)-log2(repmat(1/numel(p_m),1,numel(p_m)))));
%        ent_motifs(i) = -nansum(p_m.*log(p_m)); %shannon entropy
%        ent_motifs(i) = nansum(pdist(p_m','euclidean').^2);
       
       %compute the entropy of the transition matrix between motifs
       p_t = zeros(N+1,N+1);
       for k=1:numel(temp_idx)-1 %tally the transitions
           p_t(temp_idx(k),temp_idx(k+1)) = p_t(temp_idx(k),temp_idx(k+1)) + 1;
       end       
       p_t=p_t.*~eye(size(p_t)); %remove self transitions       
       num_transitions(i) = sum(p_t(:)); %save off the number of transitions
       p_t = (p_t(:)+eps)/sum(p_t(:)); %vectorize and change to percent
       p_t(p_t<=eps)=[]; %remove non-existent transitions
%        num_transitions(i) = numel(p_t); %number of unique transitions
       ent_transitions(i) = sum(p_t' .* (log2(p_t')-log2(repmat(1/numel(p_t),1,numel(p_t))))); %KL divergence
%        ent_transitions(i) = -nansum(p_t.*log(p_t));  %shannon entropy
       
       
    end

    %FYI to gutcheck this, just plot the histogram of probabilities during
    %the lowest and highest entropy. They should be roughly more uniformly
    %distributed (e.g. use more motifs more consistently across motifs
    %in high entropy. 
    
    
end



