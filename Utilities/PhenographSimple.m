function [ovr_comm, ovr_Q] = PhenographSimple(data, varargin)

%Define options
opts.k = 50; % the number of nearby points to use
opts.DistanceMetric = 'euclidean'; %can either be mahalanobis, seuclidean or euclidean
opts.Verbose = 1;
opts.LouvainRestarts = 1; %how many restarts to do on louvain clustering algorithm (avoids local minima)

opts = ParseOptionalInputs(opts,varargin);

%% Any points in input with zero variance should be excluded
bad_dim = (nanvar(data, [], 1) <= eps);
data = data(:, ~bad_dim);

%% Compute neighbors
N = size(data, 1);
%First build k-neighbor network
if opts.Verbose, fprintf('Finding %d closest waveforms...', opts.k); end
Idx = knnsearch(data, data, 'K', opts.k, 'Distance', opts.DistanceMetric);
if opts.Verbose, fprintf('done.\n'); end

%% Calculate network graph weights

%Now apply Jacard similarity to create connectivity graph
if opts.Verbose, fprintf('Calculating weights of network graph...\n'); end
sel_N = size(data, 1);
% w = sparse(sel_N, sel_N);
temp = cell(1,sel_N);
poolobg = parpool('local',6);
parfor i = 1:sel_N
    temp{i} = sparse(1,sel_N);
    %Find those indices with overlap
    intersect_idx = sum(ismember(Idx, Idx(i, :)), 2);
    j_ind = (intersect_idx > 0);
    temp{i}(1, j_ind) = intersect_idx(j_ind)./(2*opts.k - intersect_idx(j_ind));
    if opts.Verbose & mod(i, round(sel_N/10)) == 0, fprintf('\t%d%% done...\n', round(i/sel_N*100)); end
end %node loop
w = sparse(cat(1,temp{:}));
delete(poolobg);

%w = (w+w') - eye(size(w,1)).*diag(w); %symmetrize -- may not be necessary depending on community algorithm
if opts.Verbose, fprintf('done.\n'); end

%% Compute clusters

if opts.Verbose, fprintf('Clustering...'); end
[ovr_comm, ovr_Q] = LouvainCommunity(w, 'NumRandomStarts', opts.LouvainRestarts);
ovr_Q = ovr_Q(1); %select only the best Q
if opts.Verbose, fprintf('done.\n'); end

if opts.Verbose
    for i = 1:size(ovr_comm, 2)
        fprintf('Level %d: Found %d communities...\n', i, length(unique(ovr_comm(:, size(ovr_comm, 2) - i + 1))));
    end
end

end
