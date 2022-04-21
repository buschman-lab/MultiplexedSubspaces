function neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth,probe_num)
%THERE IS GOING OT BE A MISSMATCH IF you try to do a probes in wrong order so you
%need to include probe_num... need to fix this
%get the anatomical locations
if nargin <3; probe_num = []; end
    
ccf_path = fileparts(EphysPath);
ccf_path = load([ccf_path filesep 'probe_ccf_xvalidated.mat']);
if isempty(probe_num)
    neu_area = arrayfun(@(n) MapAnatomicalLocation(ccf_path.st,ccf_path.probe_ccf(n),st_depth{n},1),1:numel(st_depth),'UniformOutput',0);
else
    neu_area = arrayfun(@(n) MapAnatomicalLocation(ccf_path.st,ccf_path.probe_ccf(probe_num(n)),st_depth{n},1),1:numel(st_depth),'UniformOutput',0);
end

%inverse labels to match the depth plotting (inversed)
neu_area = cellfun(@(x)  x(linspace(numel(x),1,numel(x))), neu_area,'UniformOutput',0);

%switch to acroynm
neu_area = cellfun(@(x) AreaAcryonym(x,ccf_path.st), neu_area,'UniformOutput',0);

end
