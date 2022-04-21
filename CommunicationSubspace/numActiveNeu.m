function [n,n_source] = numActiveNeu(data)
%get the number of active neurons per motif per brain area

paired_areas = data(1).paired_areas;
area_label = data(1).area_label;

%get the number of neurons per area (since this can vary by motif)
n = cat(2,data(:).rrr_B);
n = cellfun(@(x) size(x,2), n);
n_source = n;
n = arrayfun(@(nn) n(paired_areas(:,2)==nn,:), 1:numel(area_label),'UniformOutput',0);
n = cellfun(@(x) x(1,:),n,'UniformOutput',0);
n = cat(1,n{:});

end %function end