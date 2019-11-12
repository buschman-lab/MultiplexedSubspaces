function [onset, offset] = ManualThreshold(trace)
% Showing the timing trace and manually select the threshold

figure; hold on; plot(trace);
threshold = inputdlg('Threshold value','Choose Timing Signal Onset Threshold'); 
threshold = str2num(threshold{1});
close

onset = find(trace>threshold,1,'first');
offset = find(trace>threshold,1,'last');

end