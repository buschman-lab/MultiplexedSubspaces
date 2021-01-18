function [dff,avgproj] = makeDFF(stack, opts, type, w)
%This function takes the imput raw stack and makes a dff along the third dimension. 
if nargin <3; type = 'dff'; end
if nargin <4; w = opts.fps*opts.method_window; end

%Choose type of averaging
if strcmp(opts.method,'mean')
    avgproj  = nanmean(stack,3);
elseif strcmp(opts.method,'median')
    avgproj  = nanmedian(stack,3);
elseif strcmp(opts.method,'mode')
    avgproj  = mode(stack,3);
elseif strcmp(opts.method,'movingavg') %30 second moving average window (since opts.fps is multiplexed
    avgproj = movmean(stack,w,3,'Endpoints','shrink');
else
    error('unknown dff method');
end

%calculate dff
dff = zeros(size(stack,1),size(stack,2),size(stack,3));
for i = 1:size(stack,3)       
    if strcmp(opts.method,'movingavg')        
        if strcmp(type,'fractional')%get fractional 
            dff(:,:,i) = ((double(stack(:,:,i)))./avgproj(:,:,i));
        else %get dff
            dff(:,:,i) = ((double(stack(:,:,i))-avgproj(:,:,i))./avgproj(:,:,i))*100;
        end
    else
        if strcmp(type,'fractional')%get fractional 
            dff(:,:,i) = (double(stack(:,:,i))./avgproj);
        else %dff
            dff(:,:,i) = ((double(stack(:,:,i))-avgproj)./avgproj)*100;
        end
    end
end

%linear detrend 
if opts.detrend
    dff = detrendNaN3(dff);
end
    
end



