function [x_alligned, onset] = AllignEphysToImaging(x,cam,frames)
%Camden MacDowell - timeless
%x and cam should be N x T matrices/vectors, frames is a start and stop
%frame index
if nargin <3; frames = []; end    

gp = general_params; 


%threshold the camera signal and times of each frame
if size(x,2)~= numel(cam) %if already LFp
    cam = interp1((1:1:numel(cam)), cam, (1:numel(cam)/size(x,2):numel(cam)), 'linear');
end

%threshold
c_bin = diff(cam > gp.camera_threshold);
onset = find(c_bin==1);  

if ~isempty(frames) %currently just locking to th eonset of each frame    
   idx = [onset(frames(1)),onset(frames(2))];
else %use entire time
   idx = [onset(1),onset(end)];
end

x_alligned = x(:,idx(1):idx(2));



end



