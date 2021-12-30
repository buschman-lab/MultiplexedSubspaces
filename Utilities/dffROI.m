function [trace,opts] = dffROI(dff,nanpxs,coords,varargin)
%Camden MacDowell - timeless
%returns a roi x time matrix of average values within a dff ROI 
%Can also be used for motif reconstructions, etc.
%coords is cell array of 2x1 coordinates of each roi

%mutable options
opts.radius = 2; %pixel radius around probe tip    
opts.offset = {[2,0],[0,-1],[1,-1],[0 0]}; %moves center of DFF circle for probes at angles. [x,y]

opts = ParseOptionalInputs(opts,varargin);

%get dimesions of the full data
x = sqrt(size(dff,2)+numel(nanpxs));   
z = size(dff,1);

%parse the radius around the probes
trace = NaN(z,numel(coords));
for cur_roi = 1:numel(coords)
    %create mask of radius around probe tip
    temp = coords{cur_roi};
    mask = zeros(x,x);  
    temp(1,1)=temp(1,1)+opts.offset{cur_roi}(1);
    temp(1,2)=temp(1,2)+opts.offset{cur_roi}(2);
    mask(temp(1,2)-opts.radius:temp(1,2)+opts.radius,temp(1,1)-opts.radius:temp(1,1)+opts.radius)=1;
    %flatten
    mask = mask(:);
    mask(nanpxs)=[];        

    trace(:,cur_roi) = nanmean(dff(:,mask==1),2);    
end    

end %function end