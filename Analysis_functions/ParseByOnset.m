function [trig_dff,trig_st] = ParseByOnset(dff,st,motif_onset,win,cur_motif)
%Camden - timeless
%takes onsets and returns the triggered onsets for all probes in dff and st

%dff input is dur x n_probe
%pad the dff probe so that you can grab edges
padsize = sum(abs(win));
if ~isempty(dff)
    dff = cat(1,NaN(padsize,size(dff,2)), dff, NaN(padsize,size(dff,2)));
    trig_dff = arrayfun(@(x) dff(x+win(1)+padsize:x+win(2)+padsize,:)',motif_onset{cur_motif},'UniformOutput',0);
    trig_dff = cat(3,trig_dff{:}); %loc x dur x onset   
else
    trig_dff = [];
end

%st input is cell array of probe. with each = dur x num units
if ~isempty(st)
    trig_st = cell(1,numel(st)); %keep as cell since diff num neurons per probe
    for i = 1:numel(st) %probe loop
       temp = st{i};
       temp = cat(1,NaN(padsize,size(temp,2)), temp, NaN(padsize,size(temp,2)));
       trig_st{i} = arrayfun(@(x) temp(x+win(1)+padsize:x+win(2)+padsize,:)',motif_onset{cur_motif},'UniformOutput',0);
       trig_st{i} = cat(3,trig_st{i}{:}); %neuron x dur x onset
    end %probe loop
else
    trig_st = [];
end

end %function end

