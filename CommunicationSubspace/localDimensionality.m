function [area_label,d_opt,d,n] = localDimensionality(data)
%Camden - timeless
%gets the dimensionality for each brain region (from factor analysis)
%during each motifs activity (rows of structure "data")
%returns the area_label, the optimal 
%number of local dim, the dim needed to capture all shared local variance
%and the number of neurons in that

d_opt = cell2mat(cat(2,data(:).qOpt_target));
%make sure none are over 30 (the limit)
if sum(d_opt==30)>0
    warning('some locations dimensionality is limited by q')
end

paired_areas = data(1).paired_areas;
area_label = data(1).area_label;

d_opt = arrayfun(@(nn) d_opt(paired_areas(:,2)==nn,:), 1:numel(area_label),'UniformOutput',0);

%gutcheck: d should be same across repeat iterations (if not, then this
%means the rng wasn't reset) ... check first motif
if sum(cellfun(@(x) sum(diff(x(:,1))), d_opt))>0
    warning('dim not consistent, check the rng seed');
end

d_opt = cellfun(@(x) x(1,:),d_opt,'UniformOutput',0);
d_opt = cat(1,d_opt{:});

%full dimensionality
d = cat(2,data(:).cvl_fa_target);
d = cellfun(@(x) find(x<0.001,1,'first'), d,'UniformOutput',1);
d = arrayfun(@(nn) d(paired_areas(:,2)==nn,:), 1:numel(area_label),'UniformOutput',0);
d = cellfun(@(x) x(1,:),d,'UniformOutput',0);
d = cat(1,d{:});

n = numActiveNeu(data); 


    
end %function end