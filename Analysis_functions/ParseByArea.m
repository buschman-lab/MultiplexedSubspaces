function [area_val, area_label] = ParseByArea(x,neu_area,type)
%Camden MacDowell - timeless
%takes any input vector of equal size to number of units and parses by the
%unique anatomical areas in neu_area

switch type
    case 'detail'
        temp = cat(1,neu_area(:).detailed_label);
    case 'parent'
        temp = cat(1,neu_area(:).parent_label);
    case 'general'
        temp = cat(1,neu_area(:).parent_label);            
        temp(ismember(temp,{'RSPagl','RSPd'}))={'RSP'}; %group RSP
        temp(ismember(temp,{'VISp','VISpm','VISa','VISam'}))={'VIS'}; %group VIS
        temp(ismember(temp,{'SSp','SSp-un','SSp-n','SSp-m','SSp-tr','SSp-ul','SSp-ll'}))={'SS'}; %group Sensory (keep bfd separate)
        temp(ismember(temp,{'MED','EPI','ILM','LAT','mfbse','MTN'}))={'THAL'}; %group Thalamo (keep bfd separate)
        temp(ismember(temp,{'CA','DG','ILM','LAT'}))={'HIPP'}; %group Hippo (keep bfd separate)
        temp(ismember(temp,{'ACAd','PL','ILA','ORBm'}))={'PRE'}; %group Front (keep MOs separate)
    otherwise
        error('unknown level of anatomical detail');
end

%split by unique areas
area_label = unique(temp);
area_val = cellfun(@(y) x(strcmp(temp,y),:,:,:), area_label,'UniformOutput',0); %supports up to 4D


end %function