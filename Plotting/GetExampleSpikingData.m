function [area_val,area_label] = GetExampleSpikingData(cur_rec,fast_load)
%Camden MacDowell - timeless
if nargin <2; fast_load=1; end
if fast_load==0
    [~,~,~,EphysPath,~] = LoadDataDirectories(cur_rec);
    [st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe'); 
else %preprocessed
    [rec_name,~,~,EphysPath] = LoadDataDirectories(cur_rec);
    [fn,~] = GrabFiles([rec_name '\w*RRR_muaflag1_GROUPEDmotif1.mat'],0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\analysisplayground\CCA'}); %first motif has all that you need
    temp = cellfun(@(x) load(x,'st_mat','st_depth'),fn);  
    st_depth = temp.st_depth;
    st_mat = temp.st_mat; 
end
%split into areas
st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
[area_val, area_label] = ParseByArea(cat(2,st_norm{:})',cat(2,neu_area{:}),'general');
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

end