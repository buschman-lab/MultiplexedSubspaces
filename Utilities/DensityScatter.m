function [fh_scat, fh_density] = DensityScatter(x,y,varargin)
% Bin the data:
opts.n_pts = 101;
opts.kernel = 1;
opts.sigma = 2.5;
opts = ParseOptionalInputs(opts,varargin);

    
pts = linspace(round(min(cat(1,x(:),y(:))),2), round(max(cat(1,x(:),y(:))),2), opts.n_pts);
N = histcounts2(y(:), x(:), pts, pts);

% Create Gaussian filter matrix:
[xG, yG] = meshgrid(-opts.kernel:opts.kernel);
g = exp(-xG.^2./(2.*opts.sigma.^2)-yG.^2./(2.*opts.sigma.^2));
g = g./sum(g(:));

% Plot scattered data (for comparison):
fh_scat = figure; hold on; 
scatter(x, y, 'r.');
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));

% Plot heatmap:
fh_density = figure; hold on; 
imagesc(pts, pts, conv2(N, g, 'same'),[0 2]);
axis equal;
colorbar
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');

end