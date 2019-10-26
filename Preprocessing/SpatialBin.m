function imds = SpatialBin(imstack,dsFactor,vmask,IgnoreNan)

% spatially bins image by a factor of dsFactor 
%LucasPinto 2017
if nargin <3
    vmask = [];
end
if nargin <4
    IgnoreNan = 1;
end

% size
[nX,nY,nZ]  = size(imstack);
imds        = zeros(ceil(nX/dsFactor),ceil(nY/dsFactor),nZ);
xbins       = unique([1:dsFactor:nX nX+1]);
ybins       = unique([1:dsFactor:nY nY+1]);

%check if input is gpuArray and make output agree
if isa(imstack,'gpuArray')
    imds = gpuArray(imds);
end
    

% downsample and ignore vasculature pixels
for zz = 1:nZ
  im        = imstack(:,:,zz);
  if ~isempty(vmask)
    im(vmask==0) = nan;
  end
  for xx = 1:length(xbins)-1
    for yy = 1:length(ybins)-1
      vals              = im(xbins(xx):xbins(xx+1)-1,ybins(yy):ybins(yy+1)-1);
      if IgnoreNan
        imds(xx,yy,zz)  = nansum(vals(:));
      else
        imds(xx,yy,zz)  = sum(vals(:));
      end
    end
  end
end