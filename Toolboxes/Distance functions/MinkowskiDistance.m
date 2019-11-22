function err = MinkowskiDistance(resp, targ, varargin),

if length(varargin) >= 1,
    p = varargin{1};
else
    p = 3;
end

err = sum((resp(:) - targ(:)).^p).^(1/p);