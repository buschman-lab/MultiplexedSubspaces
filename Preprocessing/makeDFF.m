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
elseif strcmp(opts.method,'movingavg') 
    avgproj = movmean(stack,w,3,'Endpoints','shrink');    
elseif strcmp(opts.method,'zscore') %divide by moving standard deviation
    avgproj = movmean(stack,w,3,'Endpoints','shrink');
    stdproj = movstd(stack,w,0,3,'Endpoints','shrink');
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
%             dff(:,:,i) = ((double(stack(:,:,i))-avgproj(:,:,i)));
        end
    elseif strcmp(opts.method,'zscore')
        if strcmp(type,'fractional')%get fractional 
            error('fractional dfs not support when zscoring')
        else %get dff
            dff(:,:,i) = ((double(stack(:,:,i))-avgproj(:,:,i))./stdproj(:,:,i));
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



