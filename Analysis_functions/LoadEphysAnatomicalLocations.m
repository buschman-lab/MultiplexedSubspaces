function neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth)

%get the anatomical locations
ccf_path = fileparts(EphysPath);
ccf_path = load([ccf_path filesep 'probe_ccf_xvalidated.mat']);
neu_area = arrayfun(@(n) MapAnatomicalLocation(ccf_path.st,ccf_path.probe_ccf(n),st_depth{n},1),1:numel(st_depth),'UniformOutput',0);
%inverse labels to match the depth plotting (inversed)
neu_area = cellfun(@(x)  x(linspace(numel(x),1,numel(x))), neu_area,'UniformOutput',0);
%switch to acroynm
neu_area = cellfun(@(x) AreaAcryonym(x,ccf_path.st), neu_area,'UniformOutput',0);

end
