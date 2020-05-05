function [dff_norm_lin_array, nanpxs_array, stats] = linearizeDff(dffarray,varargin)
%input dff is cell array, output dff_norm_lin is also cell array same with
%nanpxs

%Filter, spatially bins, smooths, thresholds, nans and linearizes dff for
%prep for seqNMF

opts.sm_kern = []; % [1 1 1; 1 1 1; 0 0 0];
opts.spatialbin = 2;
opts.filtband = [0.1 4];
opts.ftype = 'bandpass';
opts.nSTD = 2;
opts.Verbose = 1; 
opts.fps = 13; %original data frame rate


%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin)
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end

%%
%preallocate 
dff_norm_lin_array  = cell(1,size(dffarray,2));
nanpxs_array = cell(1,size(dffarray,2));

%loop through each dff; 
for cur_rec = 1:size(dffarray,2)
    if opts.Verbose; fprintf('\tProcessing rec %d of %d...\n',cur_rec, size(dffarray,2)); end     
    dff = dffarray{cur_rec};

    %Optional Smooth
    if ~isempty(opts.sm_kern)
%         dff = imgaussfilt3(dff,opts.sm_kern);
        dff = SpatialGaussian(dff);
%         opts.sm_kern = opts.sm_kern./sum(opts.sm_kern(:));
% 
%         for T = 1:size(dff, 3)    
%             dff(:, :, T) = nanconv(dff(:,:,T), opts.sm_kern, 'same','nanout',true);
%         end
    end
    
    %Optional Filter
    if ~isempty(opts.filtband)
        dff = filterstack(dff,opts.fps, opts.filtband, opts.ftype, 1, 0);
    end

    %Find bursts by thresholding 
    [dff, nanpxs] = conditionDffMat(dff);
    stats = struct();
    for px = 1:size(dff,2)
        temp = dff(:,px)';
%         temp(temp<=(opts.nSTD*std(temp)))=0; %Previous
%         version as of 6/10/2019
        stats.perc_points = sum(temp<=(mean(temp)+(opts.nSTD*std(temp))))/numel(temp);
        stats.thresh_val = mean(temp)+(opts.nSTD*std(temp));
        temp(temp<=(mean(temp)+(opts.nSTD*std(temp))))=0; 
        dff(:,px) = temp';
    end
    dff = conditionDffMat(dff, nanpxs);
       
    %Spatial bin
    if ~isempty(opts.spatialbin)
        dff = SpatialBin(dff,opts.spatialbin);
    end
    
    %Find bad dimensions (e.g. those with zero variance)
    dff = reshape(dff,[size(dff,1)*size(dff,2),size(dff,3)]);
    bad_dim = (nanvar(dff, [], 2) <= eps);
    dff(bad_dim,:) = NaN;
    dff = reshape(dff,[sqrt(size(dff,1)),sqrt(size(dff,1)),size(dff,2)]);
    
    %linearize and normalize
    [dff, nanpxs] = conditionDffMat(dff);
    min_dff = min(dff(:));
    max_dff = max(dff(:));
%     pct_dff = prctile(dff,95,'all');

    %OPTION 1
%     norm_lin_dff = [];
%     for px = 1:size(dff,2)
%         norm_lin_dff(:,px) = dff(:,px)/(pct_dff+max(dff(:,px)));
%     end

    %OPTION 2
    norm_lin_dff =  (dff-min_dff)./(max_dff-min_dff);
 
    
    %Store
    dff_norm_lin_array{cur_rec} = norm_lin_dff';
    nanpxs_array{cur_rec} = nanpxs;

end %dff cur_rec loop


























