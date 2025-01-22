function Figure5(ovr_data, clust_list, varargin),

% Process optional inputs
opts.ClustVar = 'raw';
opts.DistFun = 'correlation';
opts.ColorBarLabel_Offset = 0;
opts.ColorBarLabel_Height = 0.025;
opts.SortClusters = 0; %do we need to sort clusters by size?

for i = 1:2:length(varargin),
    opts.(varargin{i}) = varargin{i+1};
end

%Number of clusters
k = size(clust_list, 1);

% Create color map for clusters
%clust_colmap = lines(k);
clust_colmap = turbo(k);

% % Load some plotting functions
% if ispc,
%     addpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis\Plotting\');
% else
%     addpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/Plotting');    
% end
% Create color map for areas
fp = fig_params_cortdynamics;

% Make sure our cluster list is the right shape
for i = 1:length(clust_list),
    if size(clust_list{i}, 1) == 1,
        clust_list{i} = reshape(clust_list{i}, [size(clust_list{i}, 2) 1]);
    end
end

% Parse the distance function
if strcmpi(opts.DistFun, 'AbsCorrelation'),
    dist_metric = 'correlation';
elseif strcmpi(opts.DistFun, 'AbsCosine'),
    dist_metric = 'cosine';
else
    dist_metric = opts.DistFun;
end

%% Extract out the correlation map

if strcmpi(opts.ClustVar, 'z_raw'),
    % fisherZ transformed of the raw correlation (rho)
    orig_rho_vect = VectorizeImage(fisherZ(cat(3, ovr_data.rho)));
elseif strcmpi(opts.ClustVar, 'rho'),
    % thresholded rho
    orig_rho_vect = VectorizeImage(cat(3, ovr_data.rho_thresh));
elseif strcmpi(opts.ClustVar, 'raw'),
    % raw correlation
    orig_rho_vect = VectorizeImage(cat(3, ovr_data.rho));
elseif strcmpi(opts.ClustVar, 'prct_rho'),
    % percentile of the observed rho in the permuted distribution
    orig_rho_vect = VectorizeImage(cat(3, ovr_data.prct_rho));
elseif strcmpi(opts.ClustVar, 'zrho'),
    % z-score of the observed correlation in the permuted distribution
    orig_rho_vect = VectorizeImage(cat(3, ovr_data.z_rho));
else
    error('Bad cluster variable name.');
end

%% Plot the clusters

% Sort clusters by their size
if opts.SortClusters,
    [clust_size, sort_ind] = sort(cellfun(@(x) size(x, 1), clust_list), 1, 'descend');
    sorted_clust_list = clust_list(sort_ind);
else
    sorted_clust_list = clust_list;
    clust_size = cellfun(@(x) length(x), sorted_clust_list);
end


%% Plot an image of the similarity across all images

% Grab the vectors, sorted by clusters (which are sorted by their size)
sorted_rho_vect = orig_rho_vect(cat(1, sorted_clust_list{:}), :);

% Create distance matrix
sorted_sim_mat = pdist(sorted_rho_vect, dist_metric);
sorted_sim_mat = squareform(sorted_sim_mat);

% Normalize and then convert into a similarity matrix (pdist creates a
% dissimilarity matrix)
if strcmpi(opts.DistFun, 'cityblock'),
    sorted_sim_mat = 1 - sorted_sim_mat/size(sorted_rho_vect, 2);
elseif strcmpi(opts.DistFun, 'euclidan') || strcmpi(opts.DistFun, 'seuclidean') || strcmpi(opts.DistFun, 'mahalanobis'),
    sorted_sim_mat = -sorted_sim_mat;
elseif strcmpi(opts.DistFun, 'AbsCorrelation') || strcmpi(opts.DistFun, 'AbsCosine'),
    sorted_sim_mat = abs(1 - sorted_sim_mat);
else
    sorted_sim_mat = 1 - sorted_sim_mat;
end

% Create figure with the similarity between clusters
clust_sim_fig = figure;
imagesc(sorted_sim_mat);
hold on;
axis square;
axis off;
v = axis;
title('Similarity of Subspace Networks','Fontweight','normal')

% Create rectangles on the outside to indicate the clusters
rect_edges = [0; cumsum(clust_size)];
for i = 1:k,
    % Along bottom
    h = rectangle('Position', [(rect_edges(i)+0.5) (1+opts.ColorBarLabel_Offset)*rect_edges(end) clust_size(i) opts.ColorBarLabel_Height*rect_edges(end)], 'FaceColor', clust_colmap(i, :), 'EdgeColor', 'none');
    set(h, 'Clipping', 'off');
    text((rect_edges(i)+clust_size(i)/2+0.5), (1+opts.ColorBarLabel_Offset + opts.ColorBarLabel_Height/2)*rect_edges(end), sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
    % Along left side
    h = rectangle('Position', [-(opts.ColorBarLabel_Offset+opts.ColorBarLabel_Height)*rect_edges(end) (rect_edges(i)+0.5) opts.ColorBarLabel_Height*rect_edges(end) clust_size(i)], 'FaceColor', clust_colmap(i, :), 'EdgeColor', 'none');
    set(h, 'Clipping', 'off');
    text(-(opts.ColorBarLabel_Offset+opts.ColorBarLabel_Height/2)*rect_edges(end), (rect_edges(i)+clust_size(i)/2+0.5), sprintf('%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
    % Draw white box
    plot([rect_edges(i) rect_edges(i) rect_edges(i+1) rect_edges(i+1) rect_edges(i)], [rect_edges(i) rect_edges(i+1) rect_edges(i+1) rect_edges(i) rect_edges(i)], 'w-', 'LineWidth', 2);
end


end


function [img_vect, img_sz, good_ind] = VectorizeImage(img, varargin),

in_good_ind = [];
if nargin > 1,
    in_good_ind = varargin{1};
end

% Extract variables
img_sz = [size(img, 1) size(img, 2)];
N = size(img, 3);

% Loop here rather than trying to do some tricky concatenation to avoid a
% mistake in indexing/concatenating the wrong thing
img_vect = NaN*ones(prod(img_sz), N);
for i = 1:N,
    img_vect(:, i) = reshape(img(:, :, i), [prod(img_sz) 1]);
end

% Find the points that have no NaNs anywhere
good_ind = ~any(isnan(img_vect), 2);
if ~isempty(in_good_ind),
    good_ind = good_ind & in_good_ind;
end

% Remove the bad indices
img_vect = img_vect(good_ind, :)';

end
