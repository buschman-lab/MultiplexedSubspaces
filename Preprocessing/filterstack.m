function fStack = filterstack(dff,fs, filtband,type, verbose, addDCFlag)

%Can accept both 2D (time x pixel matrix) and 3D (dff stack); 
%Lucas pinto - modified by Camden MacDowell 2018

if nargin < 2; fs = 13; end
if nargin < 3; filtband = [1 3]; end
if nargin < 4; type = 'bandpass'; end
if nargin < 5; verbose = 1; end
if nargin <6; addDCFlag = 0; end

switch type
  case 'bandpass'
    [B,A] = butter(10,[filtband(1)/(fs/2) filtband(2)/(fs/2)],'bandpass');
  case 'highpass'
    [B,A] = butter(10,filtband(1)/(fs/2),'high');
  case 'lowpass'
    [B,A] = butter(10,filtband(2)/(fs/2),'low');
end

[nX,nY,nZ] = size(dff);

if nZ > 1
  fStack    = zeros(nX,nY,nZ);
  for xx = 1:nX
    for yy = 1:nY
        if addDCFlag
            fStack(xx,yy,:) = filtfilt(B,A,squeeze(dff(xx,yy,:)))+nanmean(dff(xx,yy,:)); % add back dc
        else
            fStack(xx,yy,:) = filtfilt(B,A,squeeze(dff(xx,yy,:)));%+nanmean(dff(xx,yy,:)); % add back dc
        end
    end
    if verbose        
        if mod(xx,round(0.1*nX)) ==0
         fprintf('\t%.2g%% done filtering...\n', xx./nX*100);
        end
    end
  end
else
  fStack    = zeros(nX,nY);
  for yy = 1:nY
      if addDCFlag
        fStack(:,yy) = filtfilt(B,A,dff(:,yy))+nanmean(dff(:,yy)); % add back dc
      else
        fStack(:,yy) = filtfilt(B,A,dff(:,yy));%+nanmean(dff(:,yy)); % add back dc
      end
      if verbose        
        if mod(yy,round(0.1*nY)) ==0
         fprintf('\t%g%% done filtering...\n', yy./nY*100);
        end
      end
  end
end


