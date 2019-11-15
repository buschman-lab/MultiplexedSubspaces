function [onset, offset] = ManualThreshold(trace)
% Showing the timing trace and manually select the threshold

figure; hold on; plot(trace);
x = []; 
while true
    loop = input('Continue? (1/0) : '); %1 to continue 0 to break
    if loop ==1
        threshold = input('Threshold Value? : ');
        onset = find(trace>threshold,1,'first');
        offset = find(trace>threshold,1,'last');
        plot([onset,offset],[threshold,threshold],'linestyle','--','marker','*','color',[1,0,0],'markersize',20);
    else
        break
    end
end
        
close
onset = find(trace>threshold,1,'first');
offset = find(trace>threshold,1,'last');

end
