function ReAddProbeCoords()
%camden macdowell - timeless. 
%add _probe_coords.mat to the options of the full dataset
%need to be in the target dir

coordsfn = GrabFiles('_probe_coords.mat',0,{pwd});
for i = 1:numel(coordsfn)
   [~, fn] = fileparts(coordsfn{i});
   fn = [erase(fn,'_probe_coords'),'_1dff_combined.mat'];
   tempcoords = load(coordsfn{i}); 
   probe_coords = tempcoords.probe_coords;
   save(fn,'probe_coords','-append')
end

end %function end