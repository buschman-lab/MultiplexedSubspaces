function window_data=WindowedDist(data,window,dist_func)
%% Get windowed average (on each row)

if nargin <3; dist_func = @(x) std(x)/sqrt(numel(x)); end 

%for each row of data loop through
NumWindows = floor(length(data)/window);
window_data = NaN(size(data,1),NumWindows);
for cur_row = 1:size(data,1)
    X = data(cur_row,:);
    %Split vector into equal sized parts
    X = X(1:(NumWindows*window)); %Remove trailing timepoints that don't fill a whole window 
    Xsplit=reshape(X,[window,NumWindows]);

    %Get average for each window
    x_window=zeros(1,size(Xsplit,2));
    for cur_W = 1:size(Xsplit,2)
        x_window(cur_W) = dist_func(Xsplit(:,cur_W));
    end %avg loop end
    window_data(cur_row,:) = x_window;
end %function end