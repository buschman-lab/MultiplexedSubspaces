function hex_h = PlotHexagonGrid(ax_h, pts, hexagons, varargin),

%Plots a grid of hexagons on figure. 

% Parse inputs
opts.PlotZ = NaN; % if a number, then plots using plot3 at that z height. Useful for plotting over images.
opts.EdgeColor = 'black'; % color of the edge
opts.EdgeStyle = '-'; % transparency of the edge
opts.EdgeWidth = 0.5;
opts.FaceColor = 'none'; % color of the face
opts.FaceAlpha = 1; % transparency of the face

% Parse optional inputs
while length(varargin) >= 2,
    opts.(varargin{1}) = varargin{2};
    varargin = varargin(3:end);
end

if isempty(ax_h),
    fig_h = figure;
    ax_h = gca;
end

% Add a z dimension?
if ~isnan(opts.PlotZ),
    pts = cat(2, pts, opts.PlotZ*ones(size(pts, 1), 1));
end

% Create the patches
hex_h = patch('Faces', hexagons, 'Vertices', pts, 'FaceColor', opts.FaceColor, 'FaceAlpha', opts.FaceAlpha, 'EdgeColor', opts.EdgeColor, 'LineStyle', opts.EdgeStyle, 'LineWidth', opts.EdgeWidth);