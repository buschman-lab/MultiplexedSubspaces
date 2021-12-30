function [onsets,onsets_mask,thresh_level,h_wf,rho,rho_null]=MotifOnset(w,h,data,rho_lvl)
%Camden MacDowell - timeless
%finds when a motif is active using adaptive thresholding
%code works, but syntax needs to be cleaned up (lots of 'temp' variables)

if nargin <4; rho_lvl = 0.2; end %minimum rho_lvl a motif has to achieve to be considered active

[p,n] = size(data); 
[~,m,L] = size(w);

%loop through each motif
rho=NaN(m,9);
rho_null=NaN(m,9);
h_wf = cell(1,m);
thresh_level=NaN(1,m);
onsets = cell(1,m);
onsets_mask = ones(1,m); %mastk to ignore motifs that are not well represented (i.e. bnelow rho_lvl)
for cur_m = 1:m
    
    %multiple thresholdlevels
    thresh = linspace(mean(h(cur_m,:)),max(h(cur_m,:)),10);    
    for t = 1:numel(thresh)-1 %threshold loop
        %threshold
        temp = h(cur_m,:);
        temp(temp<=thresh(t))=0;
        temp(temp>0)=1; 
        
        %get positive crossings
        cur_onset = find(diff(temp)==1); 
                
        %get motif-triggered raw activity
        chunk = NaN(p,L,numel(cur_onset));
        h_wf_temp = NaN(9,numel(cur_onset));
        for i = 1:numel(cur_onset)
            %for this estimate, ignore any truncated motif occurances
            if cur_onset(i)<n-L
                chunk(:,:,i) = data(:,cur_onset(i):cur_onset(i)+L-1);
            end
            
            %get the waveform of H around occurance
            if cur_onset(i)>5
                h_wf_temp(:,i) = h(cur_m,cur_onset(i)-4:cur_onset(i)+4); 
            end
        end
        h_wf{cur_m} = nanmedian(h_wf_temp,2); %save off 
        
        %only use the 'main' active part of the motif (the ~1second)
        temp = squeeze(w(:,cur_m,:));
        [~,bad_col] = mink(nanmean(temp),size(temp,2)-15,2);
        temp(:,bad_col) = [];
        avg_chunk = nanmedian(chunk,3); 
        avg_chunk(:,bad_col) = [];
        rho(cur_m,t) = corr(avg_chunk(:),temp(:));
        rho_temp = NaN(m,1);
        for q = 1:m
            if q ~= cur_m %get correlation with motif
                temp = squeeze(w(:,q,:));
                [~,bad_col] = mink(nanmean(temp),size(temp,2)-13,2);
                temp(:,bad_col) = [];
                avg_chunk = nanmedian(chunk,3);
                avg_chunk(:,bad_col) = [];
                rho_temp(q) = corr(avg_chunk(:),temp(:));                
            end
        end
        rho_null(cur_m,t) = nanmean(rho_temp);
        
    end %threshold loop
    
    %get the optimal thresold value
    [~,thresh_level(cur_m)]=max(rho(cur_m,:));
    

    %collate onsets
    temp = h(cur_m,:);
    temp(temp<=thresh(thresh_level(cur_m)))=0;
    temp(temp>0)=1;
    onsets{cur_m} = find(diff(temp)==1); 
    if nanmax(rho(cur_m,:))<rho_lvl
        onsets_mask(cur_m) = 0;
    end
end
