function [onsets,onsets_mask,thresh_level,rho_thresh]=MotifOnset(w,h,rho_frame,rho_lvl)
%Camden MacDowell - timeless
%finds when a motif is active using adaptive thresholding
%adaptive thresholding finds the threshold level at which the
%motif-triggered raw data is maximmally correlated with the motif.
%Notes: 


if nargin <5; rho_lvl = 0.2; end %minimum rho_lvl a motif has to achieve to be considered active

[~,m,L] = size(w);
numthresh = 50;

%nan pad by the length of L to get the edges
rho_frame = [rho_frame,NaN(m,L)];

%rmv edge effect at begining
h(:,1)=0;

num_onsets = NaN(m,numthresh-2); %the number of onsets at that threshold lvl
rho = NaN(m,numthresh-2); %the correlation 
thresh_level = NaN(m,1); %the chosen threshold level
rho_thresh = NaN(m,1); %the rho at the threshold level
onsets_mask = ones(m,1); %mask any not well represented motifs
onsets = cell(m,1);
for cur_m = 1:m
    %multiple thresholdlevels
    thresh = linspace(nanmedian(h(cur_m,:)),max(h(cur_m,:)),numthresh);
    thresh = thresh(2:end);
    cur_onset = cell(1,numel(thresh)-1);
    for t = 1:numel(thresh)-1 %threshold loop
        %threshold
        temp = h(cur_m,:);
        temp(temp<thresh(t))=0;
        temp(temp>0)=1; 
        
        %get positive crossings
        cur_onset{t} = find(diff(temp)==1);        
        
        %get the trajectory of spatial correlation in the data
        rho_activity = arrayfun(@(x) fisherZ(rho_frame(cur_m,x:x+L)),cur_onset{t},'UniformOutput',0);
        
        %take the peak of the mean rho as a measure of fit 
        rho(cur_m,t) = fisherInverse(max(nanmean(cat(1,rho_activity{:}))));
            
        num_onsets(cur_m,t) = numel(cur_onset{t});
    end
    %get the optimal thresold value
    [rho_thresh(cur_m),thresh_level(cur_m)]=max(rho(cur_m,:),[],2);
    if rho_thresh(cur_m)<rho_lvl %if below rho_lvl then not well represented
        onsets_mask(cur_m) = 0;
    end
    %collect onsets 
    onsets{cur_m} = cur_onset{thresh_level(cur_m)}; 
    %replace with actual value
    thresh_level(cur_m) = thresh(thresh_level(cur_m));     
end

% figure; plot peak correlation and seperatability per threshold lvl
% figure; hold on; 
% for i = 1:m
%     subplot(5,3,i); hold on; 
%     plot(rho(i,:),'k'); 
%     yyaxis right
%     plot(num_onsets(i,:),'r'); 
% end



%         %get the seperatability between that peak and each other motif -
%         not needed, since scales closely with rho (as expected)
%         sep_all = arrayfun(@(x) rho(cur_m,t)-fisherZ(rho_frame(~ismember(1:m,cur_m),x+peak)),cur_onset,'UniformOutput',0);
%         sep_all = nanmean(cat(2,sep_all{:}),2); %from each other motif in case useful in future
%         
%         %get mean seperatability 
%         sep(cur_m,t) = nanmean(sep_all);  