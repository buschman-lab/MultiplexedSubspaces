function [dff,avgproj] = makeDFF(stack, opts, w)
%This function takes the imput raw stack and makes a dff along the third dimension. 
if nargin <3; w = opts.fps*10; end

%Choose type of averaging
if strcmp(opts.method,'mean')
    avgproj  = nanmean(stack,3);
elseif strcmp(opts.method,'median')
    avgproj  = nanmedian(stack,3);
elseif strcmp(opts.method,'mode')
    avgproj  = mode(stack,3);
elseif strcmp(opts.method,'movingavg') %10 second moving average window
    avgproj = movmean(stack,w,3,'Endpoints','shrink');
else
    error('unknown dff method');
end

%calculate dff
dff = zeros(size(stack,1),size(stack,2),size(stack,3));
for i = 1:size(stack,3)
    if strcmp(opts.method,'movingavg')
        dff(:,:,i) = ((double(stack(:,:,i))-avgproj(:,:,i))./avgproj(:,:,i))*100;
    else
        dff(:,:,i) = ((double(stack(:,:,i))-avgproj)./avgproj)*100;
    end
end

%linear detrend 
if opts.detrend
    dff = detrendNaN3(dff);
end
    
end



