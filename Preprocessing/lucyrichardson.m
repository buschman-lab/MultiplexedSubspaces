function deconv = lucyrichardson(trace,fps,rise,decay,smooth_kern)
% Shell for Lucy-Richardson decvolution using kernal the approximates the
% calcium response to neural spiking rate from chen et al., 2013 (cell). 

% the algorithm works only on strictly non negative input
%rise = 150
%decay = 550

trace = trace-min(trace);

%build the kernel
dur = round(1000/fps,0);
rise_pt = round(rise/dur,0)+1; %+1 because linespace
decay_pt = round(decay/dur,0)+1; %+1 because linespace

x = [0, 1];
y = [0, 1]; 
k = 5; %this is manually fit such that the half life is ~550ms for gcamp 6f (k=5 for 13 fps)
dy = y(2) - y(1);
dx = x(2) - x(1);

%rise kernel
f = @(x) dy.*exp(k*(x + dx)) + y(2);
x = linspace(x(1), x(2), rise_pt);
y = f(x);
%normalize to zero-->1
y_rise = (y-min(y))/(max(y)-min(y));
% plot(x,y_rise)

%decay kernel
f = @(x) dy.*exp(-k*(x + dx)) + y(2);
x = linspace(x(1), x(2), decay_pt);
y = f(x);
y_decay = (y-min(y))/(max(y)-min(y));
% plot(x,y_decay)

%combine them
conv_kernel = [zeros(1,numel(y_rise)+numel(y_decay)), y_rise, y_decay]';

deconv = deconvlucy(trace,conv_kernel);

% figure; hold on; 
% plot(trace);
% plot(deconv);
%Optional Smoothing
if ~isempty(smooth_kern)
    deconv = movmean(deconv,smooth_kern,1);
end

end



