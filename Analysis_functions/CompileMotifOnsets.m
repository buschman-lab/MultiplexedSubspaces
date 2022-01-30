function [motif_onset, chunk_names] = CompileMotifOnsets(motif_fits)
%Camden - timeless
%motif_fits is the path of all of the fits
%IMPORTANT: motif_fits should be PER RECORDING. I think the code needs
%adjustment for anything beyond that. 

%re-order by chunk - FYI chunks are loaded in sequence in 2022 grabfiles so
%not actually needed; hence commented out
% %code, but always reconfirm for future data
% chunk = zeros(numel(motif_fits),1); %assumes ordered chunks (as above, this is the case)
% for cur_fit = 1:numel(motif_fits)
%    %match to test or train data
%    [~,temp_fn]= fileparts(motif_fits{cur_fit});
%    if contains(temp_fn,'test')
%        chunk(cur_fit)=1;
%    end
% end

onsets = cell(1,numel(motif_fits));
chunk_names = cell(1,numel(motif_fits));
for cur_fit = 1:numel(motif_fits)
   fit_data = load(motif_fits{cur_fit},'w','H','stats_refit');  
   
   %remove noise motifs
%    warning('removing noise motif %d');
%    fit_data.w(:,2,:) = [];
%    fit_data.H(2,:) = [];
%    fit_data.stats_refit.rho_frame_per_motif(2,:)=[];
   
   %keep track of the names for gutchecking later
   [~,temp] = fileparts(motif_fits{cur_fit});
   chunk_names{cur_fit} = extractBetween([temp,'.'],'chunk_','.');
   
   %parse the onsets
   [onsets{cur_fit},~,~] = MotifOnset(fit_data.w,fit_data.H,fit_data.stats_refit.rho_frame_per_motif);   close; 
end %fit loop

%reorder: as of 2022 1/14/2022 it is test --> train so you want to reverse
motif_onset = onsets; 
motif_onset(1:2:end) = onsets(2:2:end); 
motif_onset(2:2:end) = onsets(1:2:end); 
temp = chunk_names; 
temp(1:2:end) = chunk_names(2:2:end); 
temp(2:2:end) = chunk_names(1:2:end);
chunk_names = temp; 

%concatenate and add offset to match to true time
offset = repmat((0:numel(chunk_names)-1)*size(fit_data.H,2),size(fit_data.H,1),1);
motif_onset = cat(2,motif_onset{:});
for i = 1:size(offset,1)
    for j = 1:size(offset,2)
        motif_onset{i,j} = motif_onset{i,j}+offset(i,j);
    end
end
motif_onset = arrayfun(@(n) [motif_onset{n,:}], 1:size(motif_onset,1),'UniformOutput',0);


end

