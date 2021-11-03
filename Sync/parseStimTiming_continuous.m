function frame_index = parseStimTiming_continuous(LogFn,opts)
% Camden MacDowell - timeless
% Gets the subsequent imaging frame after start of stimulus presentation
%% parse input
opts.StimChan = 3;
opts.FrameChan = 1; 
opts.TotalChan = 4; 
opts.StimThreshold = -0.22; 
opts.FrameThreshold = 1; 
opts.RemoveHemo = 1; 

%% Load the data
completelog = fopen(LogFn,'r');
[data,~] = fread(completelog,[opts.TotalChan+1,inf],'double');
fclose(completelog);

%% get timings 
% time in seconds
t = data(1,:); 

%negative crossings = stimstart
stimtimes(:,1) =  t(diff(data(opts.StimChan+1,:)>=opts.StimThreshold)==-1);
stimtimes(:,2) =  t(diff(data(opts.StimChan+1,:)>=opts.StimThreshold)==1);

%positive crossings = exposure start
frametimes = t(diff(data(opts.FrameChan+1,:)>=opts.FrameThreshold)==1);

%remove hemoframes
if opts.RemoveHemo
   frametimes = frametimes(1:2:end); 
end

%find the frames that start soonest after stimulus onset. 
frame_index = NaN(size(stimtimes));
jitter_mag = NaN(size(stimtimes)); %Note that there will be some small jitter here because of relatively slow imaging timescale
for i = 1:size(stimtimes,1)       
   frame_index(i,1) = find(frametimes-stimtimes(i,1)>=0,1);
   frame_index(i,2) = find(frametimes-stimtimes(i,2)>=0,1);
   jitter_mag(i) = frametimes(frame_index(i))-stimtimes(i); 
end

%Jitter shold be uniform from 0 to the imaging exposure duration
histogram(jitter_mag)

end







