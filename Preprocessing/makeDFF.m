function [dff,avgproj] = makeDFF(stack, opts)
%This function takes the imput raw stack and makes a dff. 

%Choose type of averaging
if strcmp(opts.method,'mean')
    avgproj  = nanmean(stack,3);
elseif strcmp(opts.method,'median')
    avgproj  = nanmedian(stack,3);
elseif strcmp(opts.method,'mode')
    avgproj  = mode(stack,3);
elseif strcmp(opts.method,'movingavg') %10 second moving average window
    avgproj = movmean(stack,floor(opts.fps*10),3,'Endpoints','shrink');
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



