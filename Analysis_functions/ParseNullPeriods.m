function [trig_dff,trig_st] = ParseNullPeriods(dff,st,motif_onset,win,n_active_frames)
%Camden - timeless
%identifies periods when no motifs are activity 
%returns the triggered onsets of those null periods 
%for all probes in dff and st

padsize = sum(abs(win));
motif_onset = [motif_onset{:}];

%consider the n_active_frames after each motif as the active periods
t = zeros(1,size(st{1},1)+n_active_frames);
for i = 1:numel(motif_onset)
    t(motif_onset(i):(motif_onset(i)+n_active_frames-1))=1;
end

%remove pad
t(size(st{1},1)+1:end)=[];

%find indices of consecutive non-motif activity
null_str = zeros(1,padsize+1);
idx = strfind(t,null_str);
%only keep indices that are far enough apart
%this is ugly but does the trick
good_idx = NaN(1,numel(idx));
COUNT=1;
for i = 1:numel(idx)
   if i ==1 
       good_idx(COUNT) = idx(1);
       COUNT = COUNT+1;
   else
       if idx(i)-good_idx(COUNT-1)>(padsize+1) %if far enough apart, then keep
           good_idx(COUNT)=idx(i);
           COUNT = COUNT+1;
       end
   end
end
good_idx(isnan(good_idx))=[];
   
%shift by the offset in win to avoid computing baseline from motif period
if win(1)<0
    good_idx = good_idx+abs(win(1));
end

%gutcheck that everything should equal 0 in t
temp = t';
temp = cat(1,NaN(padsize,size(temp,2)), temp, NaN(padsize,size(temp,2)));
trig_st = arrayfun(@(x) temp(x+win(1)+padsize:x+win(2)+padsize,:)',good_idx,'UniformOutput',0);
trig_st = cat(3,trig_st{:}); %neuron x dur x onset
assert(sum(trig_st(:)==1)==0,'motif period incorportated into null')
       
%dff input is dur x n_probe
%pad the dff probe so that you can grab edges
padsize = sum(abs(win));
if ~isempty(dff)
    dff = cat(1,NaN(padsize,size(dff,2)), dff, NaN(padsize,size(dff,2)));
    trig_dff = arrayfun(@(x) dff(x+win(1)+padsize:x+win(2)+padsize,:)',good_idx,'UniformOutput',0);
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
       trig_st{i} = arrayfun(@(x) temp(x+win(1)+padsize:x+win(2)+padsize,:)',good_idx,'UniformOutput',0);
       trig_st{i} = cat(3,trig_st{i}{:}); %neuron x dur x onset
    end %probe loop
else
    trig_st = [];
end

end %function end
