function [shared, not_shared] = SharedPixels(W,data)
%Nan out pixels in the test data that have zero variance in either the test or training data
assert(size(W,1)==size(data,1),'Cannot compute shared pixels on different sizes')
if size(W,3)>1 && sum(sum(nanvar(W, [], 3),2)) >=eps %this is written this way to deal with motifs that are just the repmat average frame
    data(sum(nanvar(W, [], 3),2) <= eps) = NaN;
else
    data(nanvar(squeeze(W(:,:,1)),[],2) <= eps) = NaN;
end
data(nanvar(data, [], 2) <= eps) = NaN;

%Remove all NaN pixel from both testing and training data
not_shared = isnan(data(:,1));
shared = ~isnan(data(:,1));
end