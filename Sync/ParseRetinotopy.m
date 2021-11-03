function frame_index = ParseRetinotopy(dff,frame_index,StimOrder_Fn,varargin)
% Camden MacDowell - timeless

% frame_index = parseStimTiming(LogFn);

opts.baselineframes = 15; 
opts.nsamples = max(frame_index(:,2)-frame_index(:,1)); %take the longest stimuli

%compute stimulus response
dff_trial = NaN(size(dff,1),size(dff,2),opts.nsamples+opts.baselineframes,size(frame_index,1));
for i = 1:size(frame_index,1)
    temp = dff(:,:,frame_index(i)-opts.baselineframes:frame_index(i)+opts.nsamples-1);
    %subtract baseline which is the frames before index
    dff_trial(:,:,:,i) = temp-nanmean(temp(:,:,1:opts.baselineframes),3);
end

%get average response to each cardinal diration
stim_type = load('Z:\Rodent Data\Wide Field Microscopy\ControlExperiments_WidefieldData\Retinotopy_2\stimInfo.mat','stim_type');
stim_type = stim_type.stim_type;
unique_stim_type = unique(stim_type);
dff_avg = NaN(size(dff_trial,1),size(dff_trial,2),size(dff_trial,3),numel(unique_stim_type));
for i = 1:numel(unique_stim_type)
   dff_avg(:,:,:,i) = nanmean(dff_trial(:,:,:,stim_type == unique_stim_type(i)),4); 
end


%%
for i = 1:size(dff_avg,4)
temp = dff_avg(:,:,:,i);
for j = 1:size(temp,3)
imagesc(temp(:,:,j),[-1 1]);
title(sprintf('%d frame %d',i,j-15));
colormap magma
pause(0.05);
end
end
%%

end







