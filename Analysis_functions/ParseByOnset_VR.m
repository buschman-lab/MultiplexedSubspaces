function [trig_dff,trig_st] = ParseByOnset_VR(dff,st,motif_onset,win,cur_motif,offset,num_motif_back)
%Camden - timeless
%takes onsets and returns the triggered onsets for all probes in dff and st

%dff input is dur x n_probe
%pad the dff probe so that you can grab edges
%V test to check one motif back %works because motifs at start are removed from analysis
%v_onset = motif_onset{1,cur_motif}(1+num_motif_back:end); %motif_onset{1,cur_motif}(1:end-num_motif_back)];
v_onset=motif_onset{cur_motif};
if num_motif_back>0
    %v_onset=zeros(1,length(motif_onset{cur_motif}));
    cats = repelem(1:14,[cellfun('size',motif_onset(1:14),2)]).';
    mot_list = cell2mat(motif_onset(1:14));    
    for l = 1:length(motif_onset{cur_motif})
        temp = motif_onset{cur_motif}(1,l);
        inds = temp-mot_list;
        inds(inds<=0) = inf;
        [c,idx] = sort(inds);
        prev_onset = motif_onset{1,cats(idx(num_motif_back))}(motif_onset{1,cats(idx(num_motif_back))}==temp-c(num_motif_back));
        if ~isempty(prev_onset)
            v_onset(1,l) = prev_onset;
        end
    end
end
v_onset = v_onset-offset;
v_onset(v_onset<0)=[];

padsize = sum(abs(win));
if ~isempty(dff)
    dff = cat(1,NaN(padsize,size(dff,2)), dff, NaN(padsize,size(dff,2)));
    trig_dff = arrayfun(@(x) dff(x+win(1)+padsize:x+win(2)+padsize,:)',v_onset,'UniformOutput',0);
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
       trig_st{i} = arrayfun(@(x) temp(x+win(1)+padsize:x+win(2)+padsize,:)',v_onset,'UniformOutput',0);
       trig_st{i} = cat(3,trig_st{i}{:}); %neuron x dur x onset
    end %probe loop
else
    trig_st = [];
end

end %function end

