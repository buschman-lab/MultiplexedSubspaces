function [clust_ind, ovr_Q] = Phenograph_CM(data, ts,NN,LL, varargin)

%Define options
opts.k = 4; % the number of nearby points to use
opts.DistanceMetric = 'mahalanobis'; %can either be mahalanobis, seuclidean or euclidean
opts.UseDerivative = 0; %whether to use derivative of waveform instead of raw waveform
opts.TimeDiscountTau = Inf; %discounts weights over time. Seems as if t_weight works better
opts.tWeight = 2; %relative weighting of time to waveform -- roughly how many seconds to equal average distance to overall mean
opts.Verbose = 1;
opts.SubSample = 100; %how many samples to use for initial clustering.  If 0, no subsampling used
opts.SubSampleType = 'decimate'; %If 'decimate', takes every N samples to get SubSample points; if 'random', takes random subsample
opts.LouvainRestarts = 100; %how many restarts to do on louvain clustering algorithm (avoids local minima)
opts.GeneralizationApproach = 'nearest'; %If 'nearest', takes mode unit classification of K nearest neighbors. If 'prototype', selects the nearest prototype, if 'likelihood' chooses the maximum likelihood.
opts.AutoEncSmooth = 0; %whether to use an auto-encoder to smooth inputs
opts.AutoEncReducDim = 0; %whether to use auto-encoder to reduce dimensionality
opts.AutoEnc_HiddenLayer = 8; % number of hidden units
opts.AutoEnc_L2WeightRegularization = 0.001; %reduce weights in autoencoder
opts.AutoEnc_SparsityProportion = 0.05; %average proportion of responses in hidden layer to any input, lower = sparser
opts.AutoEnc_SparsityRegularization = 1.6; %how much to weight sparsity

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end

%% Any points in input with zero variance should be excluded
% bad_dim = (nanvar(data, [], 1) <= eps);
% data = data(:, ~bad_dim);

%% Take first derivative

if opts.UseDerivative,
    data = cat(2, data, diff(data, 1, 2));
end

%% If necessary, subsample the data
N = size(data, 1);

% Get sub sample selection index
if (opts.SubSample > 0) & (N > opts.SubSample),
    sel_ind = false(N, 1);
    if strcmpi(opts.SubSampleType, 'decimate'),
        sel_ind(round([1:(N./opts.SubSample):N])) = true;
    elseif strcmpi(opts.SubSampleType, 'random'),
        sel_ind(randsample(N, opts.SubSample, 0)) = true;
    end
    %Select the chosen data points
    unsel_data = data(~sel_ind, :);
%     unsel_ts = ts(~sel_ind, :);
    data = data(sel_ind, :);
%     ts = ts(sel_ind, :);
else
    unsel_data = [];
    unsel_ts = [];
end

%% Apply auto-encoder

if opts.AutoEncSmooth || opts.AutoEncReducDim,
    %Train an autodecoder
    if opts.Verbose, fprintf('Training autoencoder...\n'); end
    autoenc = trainAutoencoder(data', opts.AutoEnc_HiddenLayer, 'ShowProgressWindow', 0, 'MaxEpochs', 2000, ...
        'L2WeightRegularization', opts.AutoEnc_L2WeightRegularization, ...
        'SparsityProportion', opts.AutoEnc_SparsityProportion, ...
        'SparsityRegularization', opts.AutoEnc_SparsityRegularization);
    
    % Smooth data?
    if opts.AutoEncSmooth,
        data = predict(autoenc, data')';
        unsel_data = predict(autoenc, unsel_data')';
    elseif opts.AutoEncReducDim,
        data = encode(autoenc, data')';
        unsel_data = encode(autoenc, unsel_data')';
    end
else
    autoenc = [];
end

%% Compute neighbors

%First build k-neighbor network
if opts.Verbose, fprintf('Finding %d closest waveforms...', opts.k); end
if opts.tWeight > 0,
    if strcmpi(opts.DistanceMetric, 'mahalanobis'),
        %Create covariance matrix, with weighting for time
        cov_mat = cov(data);
        tWeight = opts.tWeight/sum(diag(inv(cov_mat)));
        cov_mat(end+1, end+1) = tWeight;
        % Perform nearest neighbor search
        Idx = knnsearch(cat(2, data, ts), cat(2, data, ts), 'K', opts.k, 'Distance', 'mahalanobis', 'Cov', cov_mat);
    elseif strcmpi(opts.DistanceMetric, 'seuclidean'),
        %Correct t_weight for the knnsearch algorithm
        tWeight = opts.tWeight/sum(1./(nanstd(data, [], 1).^2));
        dist_scale = cat(2, nanstd(data, [], 1), sqrt(tWeight));
        % Perform nearest neighbor search
        Idx = knnsearch(cat(2, data, ts), cat(2, data, ts), 'K', opts.k, 'Distance', 'seuclidean', 'Scale', dist_scale);
    elseif strcmpi(opts.DistanceMetric, 'euclidean'),
        %Correct t_weight for the knnsearch algorithm
        tWeight = opts.tWeight/sum(size(data, 2));
        dist_scale = cat(2, ones(1, size(data, 2)), sqrt(tWeight));
        % Perform nearest neighbor search
        Idx = knnsearch(cat(2, data, ts), cat(2, data, ts), 'K', opts.k, 'Distance', 'seuclidean', 'Scale', dist_scale);
    else
        error('Unknown distance metric passed.');
    end
else
    if strcmpi(opts.DistanceMetric,'dtwDist')
        %Dynamic Time Warp Distance Function
        num = size(data,1);
        dtwDist = @(x,y)cellfun(@(y)dtw((reshape(x,[NN,LL])),y),...
            cellfun(@(y)reshape(y,[NN,LL]),(num2cell(y,2)),'UniformOutput',0));
        %Perform nearest neighbor search
        [Idx, dist] = knnsearch(data, data, 'K', opts.k, 'Distance', dtwDist); 
    elseif strcmpi(opts.DistanceMetric,'xcorrDist')
        num = size(data,1);
        xcorrDist = @(x,y)cellfun(@(y)abs(1-max(max(normxcorr2((reshape(x,[NN,LL])),y)))),...
            cellfun(@(y)reshape(y,[NN,LL]),(num2cell(y,2)),'UniformOutput',0));
        %Perform nearest neighbor search
        [Idx, dist] = knnsearch(data, data, 'K', opts.k, 'Distance', xcorrDist); 
    elseif strcmpi(opts.DistanceMetric,'corrDist')
        num = size(data,1);
        corrDist = @(x,y)cellfun(@(y)abs(1-max(max(corr2((reshape(x,[NN,LL])),y)))),...
            cellfun(@(y)reshape(y,[NN,LL]),(num2cell(y,2)),'UniformOutput',0));
        %Perform nearest neighbor search
        [Idx, dist] = knnsearch(data, data, 'K', opts.k, 'Distance', corrDist);       
    elseif strcmpi(opts.DistanceMetric,'rmseDist')
        %Pregenerate To Save Time: 
        %Generate a random sine wave as 'H'
        fs = 13; % Sampling frequency (samples per second)
        dt = 1/fs; % seconds per sample
        StopTime = 30; % seconds
        t = (0:dt:StopTime)'; % seconds
        F = 1; % Sine wave frequency (hertz)
        H = (sin(2*pi*F*t)')*0.75;
        H(H(1,:)<0)=0;
%         x = zeros(1,length(t)*5);
%         x(1,10:length(H)+9) = H(1,:);
%         H = [x x];
        %Same with the entire reconstruction
        XSynth = SynthReconstructions(data,H,NN,LL);
        
        rmseDistFun = @(x,y)rmseDist(x,y,NN,LL,H,XSynth);        
        %Perform nearest neighbor search
        [Idx, dist] = knnsearch(data, data, 'K', opts.k, 'Distance', rmseDistFun);            
    else
        %Perform nearest neighbor search
        Idx = knnsearch(data, data, 'K', opts.k, 'Distance', opts.DistanceMetric);
    end
end
if opts.Verbose, fprintf('done.\n'); end

%% Calculate network graph weights

%Now apply Jacard similarity to create connectivity graph
if opts.Verbose, fprintf('Calculating weights of network graph...\n'); end
sel_N = size(data, 1);
w = sparse(sel_N, sel_N);
for i = 1:sel_N,
    %Find those indices with overlap
    intersect_idx = sum(ismember(Idx, Idx(i, :)), 2);
    j_ind = (intersect_idx > 0);
    if ~isempty(ts),
        w(i, j_ind) = intersect_idx(j_ind)./(2*opts.k - intersect_idx(j_ind)).*exp(-abs(ts(i) - ts(j_ind))./opts.TimeDiscountTau);
    else
        w(i, j_ind) = intersect_idx(j_ind)./(2*opts.k - intersect_idx(j_ind));
    end
    if opts.Verbose & mod(i, round(sel_N/10)) == 0, fprintf('\t%d%% done...\n', round(i/sel_N*100)); end
end %node loop
%w = (w+w') - eye(size(w,1)).*diag(w); %symmetrize -- may not be necessary depending on community algorithm
if opts.Verbose, fprintf('done.\n'); end

%% Compute clusters

if opts.Verbose, fprintf('Clustering...'); end
[ovr_comm, ovr_Q] = LouvainCommunity(w, 'NumRandomStarts', opts.LouvainRestarts);
ovr_Q = ovr_Q(1); %select only the best Q
if opts.Verbose, fprintf('done.\n'); end

if opts.Verbose,
    for i = 1:size(ovr_comm, 2),
        fprintf('Level %d: Found %d communities...\n', i, length(unique(ovr_comm(:, size(ovr_comm, 2) - i + 1))));
    end
end


%% Generalize clustering to unselected data points
if (opts.SubSample > 0) %& ~isempty(unsel_ts),
    %Need to generalize to non-clustered data
    if opts.Verbose, fprintf('Generalizing clustered waveforms.\n'); end
    clust_ind = NaN*ones(N, size(ovr_comm, 2));
    for i = 1:size(ovr_comm, 2),
        if any(strcmpi(opts.GeneralizationApproach, {'likelihood', 'prototype'})),
            %Need to invert time discount back into seconds
            new_id = GeneralizeWaveformClassification(data, [], ovr_comm(:, i), unsel_data, [], ...
                NN,LL,'Approach', opts.GeneralizationApproach, 'Verbose', opts.Verbose, 'TimeDiscount', 0);
        elseif strcmpi(opts.GeneralizationApproach, 'nearest'),
            new_id = GeneralizeWaveformClassification(data, [], ovr_comm(:, i), unsel_data, [], ...
                NN,LL,'Approach', opts.GeneralizationApproach, 'Verbose', opts.Verbose, 'TimeDiscount', 0, ...
                'k', opts.k, 'DistanceMetric', opts.DistanceMetric);
        end
        clust_ind(sel_ind, i) = ovr_comm(:, i);
        clust_ind(~sel_ind, i) = new_id;
    end
else
    clust_ind = ovr_comm(:, 1);
end
