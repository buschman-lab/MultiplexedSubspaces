function [yall, yneu] = GetLocalTrialAverageStrength(area_val,area_label,type,normalizeflag,normtype)
%Camden - timeless
%get the trial average activity within each brain area during each motif

if nargin <3; type = 'mean'; end
if nargin <4; normalizeflag=0; end

%if want to rerun to get input, look at ComunicationSubspace_MotifTriggered

yneu = cell(1,numel(area_label));
for i = 1:numel(area_label)
    y = area_val{strcmp(area_label,area_label{i}),:};

    %normalize to baseline
    y = normalizeToBaseline(y,[1:2],normtype);

    %use post stimulus
    y = y(:,3:end,:);

    %get the psth
    y = nanmean(y,3);
    
    switch type
        case 'mean'
            yneu{i} = nanmean(y,2);
        case 'max'
            yneu{i} = nanmax(y,[],2);
    end
    
end %subspace identification loop

%overall average
idx = arrayfun(@(n) ones(numel(yneu{n}),1)*n, 1:numel(yneu),'UniformOutput',0);
idx = cat(1,idx{:});
temp_yneu = yneu;
if normalizeflag == 1 %zscore across areas
   yneu = cat(1,yneu{:});
   yneu = zscore(yneu);
   uidx = unique(idx);
   for j = 1:numel(uidx)
      temp_yneu{j} = yneu(idx==uidx(j)); 
   end
   yneu= temp_yneu;
elseif normalizeflag == 2 %normalize to maximum activity across areas
   yneu = cat(1,yneu{:});
   yneu = yneu/max(yneu);
   uidx = unique(idx);
   for j = 1:numel(uidx)
      temp_yneu{j} = yneu(idx==uidx(j)); 
   end
   yneu= temp_yneu;    
elseif normalizeflag == 3 %normalize to min max across areas
   yneu = cat(1,yneu{:});
   yneu = (yneu-min(yneu))/(max(yneu)-min(yneu));
   uidx = unique(idx);
   for j = 1:numel(uidx)
      temp_yneu{j} = yneu(idx==uidx(j)); 
   end
   yneu= temp_yneu;
end

yall = cellfun(@(x) nanmean(x), yneu);

end %function


























